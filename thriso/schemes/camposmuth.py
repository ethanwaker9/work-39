from ..meter import METER
from ..nizk import SigmaGAIP, PVP, poly_eval
from ..sharing import SubgroupShamir, lagrange_at_zero
from .base import Scheme


class CamposMuth(Scheme):
    name = "cm22"
    label = "CM22"
    active = True
    robust = True
    assumptions = "MT-GAIP, dGAIP"
    model = "QROM"

    def __init__(self, ga, n, t, K=16, zk_general=112, zk_special=71, transform="unruh"):
        super().__init__(ga, n, t, K, zk_general, zk_special, transform)
        self.sss = SubgroupShamir(ga.N, ga.factors, n, t)
        self.sigma = SigmaGAIP(ga, zk_general, zk_special, transform)
        self.pvp = PVP(ga, zk_general, self.sss.M)

    def keygen(self):
        METER.set_phase("keygen")
        METER.set_party(0)
        ga = self.ga
        M = self.sss.M
        self.pk = []
        self.state = {i: [] for i in range(1, self.n + 1)}
        for l in range(self.K):
            s = self.sss.secret()
            self.pk.append(ga.act(ga.E0, self.sss.exponent(s)))
            sh = self.sss.share(s)
            polys = {}
            for i in range(1, self.n + 1):
                polys[i] = [sh[i]] + [self.sss.secret() for _ in range(self.t - 1)]
            for i in range(1, self.n + 1):
                sub = {j: poly_eval(polys[j], i, M) for j in range(1, self.n + 1)}
                self.state[i].append({"share": sh[i], "poly": polys[i], "sub": sub})
        METER.send(curves=self.K)

    def secret_bits_per_party(self):
        return self.K * (self.t + self.n) * self.sss.share_bits()

    def setup(self, S, cheat=None):
        return True

    def sign(self, msg, S, cheat=None):
        METER.set_phase("sign")
        ga = self.ga
        N = ga.N
        M = self.sss.M
        scale = self.sss.scale
        S = list(S)
        T = self.T
        L = lagrange_at_zero(sorted(S), M)
        E = [ga.E0] * T
        store = {}
        for pos, i in enumerate(S):
            METER.begin_round()
            METER.set_party(i)
            rows = []
            newE = []
            for l in range(T):
                b = [self.sss.secret() for _ in range(self.t)]
                R = ga.act(ga.E0, ga.rand())
                Rp = ga.act(R, self.sss.exponent(b[0]))
                xs = {j: poly_eval(b, j, M) for j in S}
                main, pieces = self.pvp.prove(b, [(R, Rp)], xs, ("cm", i, l), scale=scale)
                applied = self.sss.exponent(b[0])
                if cheat is not None and cheat.get("party") == i and cheat.get("type") == "chain" and l == cheat.get("position", 0):
                    applied = (applied + 1) % N
                e2 = ga.act(E[l], applied)
                zk = self.sigma.prove([(R, Rp), (E[l], e2)], self.sss.exponent(b[0]), ("cmzk", i, l))
                METER.send(curves=3)
                rows.append((b, R, Rp, main, pieces, zk, E[l], e2))
                newE.append(e2)
            METER.end_round()
            METER.begin_round()
            for j in S:
                if j == i:
                    continue
                METER.set_party(j)
                for l in range(T):
                    b, R, Rp, main, pieces, zk, e1, e2 = rows[l]
                    if not self.sigma.verify([(R, Rp), (e1, e2)], zk, ("cmzk", i, l)):
                        METER.end_round()
                        return None
            METER.end_round()
            store[i] = rows
            E = newE
        chal = self.fish.challenge(E, msg)
        z = [0] * T
        METER.begin_round()
        for i in S:
            METER.set_party(i)
            for l in range(T):
                c = chal[l]
                b = store[i][l][0]
                if c != 0:
                    key = self.state[i][abs(c) - 1]
                    sg = 1 if c > 0 else -1
                    zp = [(bj - sg * L[i] * fj) % M for bj, fj in zip(b, key["poly"])]
                else:
                    zp = list(b)
                store[i][l] = store[i][l] + (zp,)
                z[l] = (z[l] + self.sss.exponent(zp[0])) % N
            METER.send(scalars=T * self.t)
        METER.end_round()
        METER.begin_round()
        for i in S:
            for j in S:
                if j == i:
                    continue
                METER.set_party(j)
                for l in range(T):
                    b, R, Rp, main, pieces, zk, e1, e2, zp = store[i][l]
                    c = chal[l]
                    val = poly_eval(zp, j, M)
                    if c != 0:
                        sg = 1 if c > 0 else -1
                        val = (val + sg * L[i] * self.state[j][abs(c) - 1]["sub"][i]) % M
                    if not self.pvp.verify_piece(j, val, main, pieces[j], [(R, Rp)], ("cm", i, l), scale=scale):
                        METER.end_round()
                        return None
                    if not self.pvp.verify_piece(0, None, main, pieces[0], [(R, Rp)], ("cm", i, l), scale=scale):
                        METER.end_round()
                        return None
        METER.end_round()
        return (chal, z)
