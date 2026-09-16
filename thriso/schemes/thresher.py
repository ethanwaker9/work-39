from ..meter import METER
from ..nizk import PVP, NITZKShamir, poly_eval
from ..sharing import SubgroupShamir
from .base import Scheme, exceptional_multipliers


class ThreshERSharK(Scheme):
    name = "thresher"
    label = "ABCP23b"
    active = True
    robust = True
    assumptions = "MT-GAIP, C_k-VPwAI"
    model = "QROM"

    def __init__(self, ga, n, t, K=16, zk_general=112, zk_special=71, transform="unruh"):
        super().__init__(ga, n, t, K, zk_general, zk_special, transform)
        self.deg = (n - 1) // 2
        self.sss = SubgroupShamir(ga.N, ga.factors, n, self.deg + 1, bound=2 * K + 1)
        self.M = self.sss.M
        self.pvp = PVP(ga, zk_general, self.M)
        small = min(p for p in ga.factors if p > n)
        self.nitzk = NITZKShamir(self.M, 128, small)
        self.mult = exceptional_multipliers(K, self.M)

    def keygen(self):
        METER.set_phase("keygen")
        METER.set_party(0)
        ga = self.ga
        M = self.M
        self.polys = {i: [self.sss.secret() for _ in range(self.deg + 1)] for i in range(1, self.n + 1)}
        self.x = sum(self.polys[i][0] for i in self.polys) % M
        self.pk = [ga.act(ga.E0, self.sss.exponent(c * self.x)) for c in self.mult]
        self.view = {j: {i: poly_eval(self.polys[i], j, M) for i in self.polys} for j in range(1, self.n + 1)}
        METER.send(curves=self.K)

    def secret_bits_per_party(self):
        return (self.deg + 1 + self.n) * self.sss.share_bits()

    def setup(self, S, cheat=None):
        return True

    def dkg_nonces(self, Q, T, cheat=None):
        ga = self.ga
        M = self.M
        polys = {}
        mains = {}
        pieces = {}
        METER.begin_round()
        for i in Q:
            METER.set_party(i)
            polys[i] = [[self.sss.secret() for _ in range(self.deg + 1)] for _ in range(T)]
            mains[i] = []
            pieces[i] = []
            for l in range(T):
                xs = {j: poly_eval(polys[i][l], j, M) for j in Q}
                mn, pc = self.nitzk.prove(polys[i][l], xs, ("vss", i, l))
                mains[i].append(mn)
                pieces[i].append(pc)
        METER.end_round()
        METER.begin_round()
        for i in Q:
            for j in Q:
                if j == i:
                    continue
                METER.set_party(j)
                for l in range(T):
                    x = poly_eval(polys[i][l], j, M)
                    if not self.nitzk.verify(j, x, mains[i][l], pieces[i][l][j], ("vss", i, l)):
                        METER.end_round()
                        return None
        METER.end_round()
        F = [ga.E0] * T
        for i in Q:
            METER.begin_round()
            METER.set_party(i)
            newF = []
            proofs = []
            for l in range(T):
                w = polys[i][l][0]
                applied = self.sss.exponent(w)
                if cheat is not None and cheat.get("party") == i and cheat.get("type") == "chain" and l == cheat.get("position", 0):
                    applied = (applied + 1) % ga.N
                f2 = ga.act(F[l], applied)
                xs = {j: poly_eval(polys[i][l], j, M) for j in Q}
                mn, pc = self.pvp.prove(polys[i][l], [(F[l], f2)], xs, ("pk", i, l), scale=self.sss.scale)
                newF.append(f2)
                proofs.append((mn, pc))
            METER.send(curves=T)
            METER.end_round()
            METER.begin_round()
            for j in Q:
                if j == i:
                    continue
                METER.set_party(j)
                for l in range(T):
                    mn, pc = proofs[l]
                    x = poly_eval(polys[i][l], j, M)
                    if not self.pvp.verify_piece(j, x, mn, pc[j], [(F[l], newF[l])], ("pk", i, l), scale=self.sss.scale):
                        METER.end_round()
                        return None
                    if not self.pvp.verify_piece(0, None, mn, pc[0], [(F[l], newF[l])], ("pk", i, l), scale=self.sss.scale):
                        METER.end_round()
                        return None
            METER.end_round()
            F = newF
        return F, polys

    def sign(self, msg, S, cheat=None):
        METER.set_phase("sign")
        ga = self.ga
        M = self.M
        Q = list(range(1, self.n + 1))
        T = self.T
        res = self.dkg_nonces(Q, T, cheat)
        if res is None:
            return None
        F, polys = res
        chal = self.fish.challenge(F, msg)
        z = [0] * T
        resp = {}
        METER.begin_round()
        for i in Q:
            METER.set_party(i)
            resp[i] = []
            for l in range(T):
                c = chal[l]
                if c != 0:
                    m = (1 if c > 0 else -1) * self.mult[abs(c) - 1]
                    r = [(bj - m * sj) % M for bj, sj in zip(polys[i][l], self.polys[i])]
                else:
                    r = list(polys[i][l])
                resp[i].append(r)
                z[l] = (z[l] + self.sss.exponent(r[0])) % ga.N
            METER.send(scalars=T * (self.deg + 1))
        METER.end_round()
        METER.begin_round()
        for i in Q:
            for j in Q:
                if j == i:
                    continue
                METER.set_party(j)
                for l in range(T):
                    c = chal[l]
                    lhs = poly_eval(resp[i][l], j, M)
                    rhs = poly_eval(polys[i][l], j, M)
                    if c != 0:
                        m = (1 if c > 0 else -1) * self.mult[abs(c) - 1]
                        rhs = (rhs - m * self.view[j][i]) % M
                    if lhs != rhs:
                        METER.end_round()
                        return None
        METER.end_round()
        return (chal, z)
