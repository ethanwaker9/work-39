from ..hashing import commit, LAMBDA_BYTES
from ..meter import METER
from ..nizk import SigmaGAIP
from ..sharing import Additive, SubgroupShamir
from .base import Scheme, exceptional_multipliers


class Sashimi(Scheme):
    name = "sashimi"
    label = "CS20"
    active = True
    robust = False
    assumptions = "MT-GAIP, dGAIP"
    model = "ROM"

    def __init__(self, ga, n, t, K=16, zk_general=112, zk_special=71, transform="unruh"):
        super().__init__(ga, n, t, K, zk_general, zk_special, transform)
        self.sigma = SigmaGAIP(ga, zk_general, zk_special, transform)
        self.sss = Additive(ga.N, t)

    def keygen(self):
        METER.set_phase("keygen")
        METER.set_party(0)
        ga = self.ga
        self.secrets = [self.sss.secret() for _ in range(self.K)]
        self.pk = [ga.act(ga.E0, s) for s in self.secrets]
        shares = [self.sss.share(s) for s in self.secrets]
        self.shares = {i: [shares[l][i] for l in range(self.K)] for i in range(1, self.t + 1)}
        METER.send(curves=self.K)

    def secret_bits_per_party(self):
        return self.K * self.sss.share_bits()

    def setup(self, S, cheat=None):
        return True

    def key_term(self, i, c, S):
        return (1 if c > 0 else -1) * self.shares[i][abs(c) - 1]

    def grp_action(self, S, T, cheat=None):
        ga = self.ga
        N = ga.N
        ctx = ("grpaction",)
        nonces = {}
        EP = {}
        pi1 = {}
        coms = {}
        METER.begin_round()
        for i in S:
            METER.set_party(i)
            nonces[i] = [ga.rand() for _ in range(T)]
            EP[i] = [ga.act(ga.E0, nonces[i][k]) for k in range(T)]
            pi1[i] = [self.sigma.prove([(ga.E0, EP[i][k])], nonces[i][k], (ctx, i, k, 1)) for k in range(T)]
            coms[i] = commit(("ep", EP[i]))
            METER.send(hashbytes=LAMBDA_BYTES)
        METER.end_round()
        METER.begin_round()
        for i in S:
            METER.set_party(i)
            METER.send(curves=T, hashbytes=LAMBDA_BYTES)
            for j in S:
                if j == i:
                    continue
                for k in range(T):
                    if not self.sigma.verify([(ga.E0, EP[j][k])], pi1[j][k], (ctx, j, k, 1)):
                        METER.end_round()
                        return None
        METER.end_round()
        E = [ga.E0] * T
        for pos, i in enumerate(S):
            METER.begin_round()
            METER.set_party(i)
            newE = []
            proofs = []
            for k in range(T):
                w = nonces[i][k]
                applied = w
                if cheat is not None and cheat.get("party") == i and cheat.get("type") == "chain" and k == cheat.get("position", 0):
                    applied = (w + 1) % N
                e2 = ga.act(E[k], applied)
                newE.append(e2)
                proofs.append(self.sigma.prove([(ga.E0, EP[i][k]), (E[k], e2)], w, (ctx, i, k, 2)))
            METER.send(curves=T)
            METER.end_round()
            METER.begin_round()
            for j in S:
                if j == i:
                    continue
                METER.set_party(j)
                for k in range(T):
                    if not self.sigma.verify([(ga.E0, EP[i][k]), (E[k], newE[k])], proofs[k], (ctx, i, k, 2)):
                        METER.end_round()
                        return None
            METER.end_round()
            E = newE
        return E, nonces

    def sign(self, msg, S, cheat=None):
        METER.set_phase("sign")
        ga = self.ga
        N = ga.N
        S = list(S)
        T = self.T
        res = self.grp_action(S, T, cheat)
        if res is None:
            return None
        E, nonces = res
        chal = self.fish.challenge(E, msg)
        z = [0] * T
        METER.begin_round()
        for i in S:
            METER.set_party(i)
            for k in range(T):
                c = chal[k]
                v = nonces[i][k]
                if c != 0:
                    v -= self.key_term(i, c, S)
                z[k] = (z[k] + v) % N
            METER.send(scalars=T)
        METER.end_round()
        return (chal, z)


class CSISharKActive(Sashimi):
    name = "csishark"
    label = "ABCP23a"
    active = True
    assumptions = "MT-GAIP, C_k-VPwAI"
    model = "QROM"

    def __init__(self, ga, n, t, K=16, zk_general=112, zk_special=71, transform="unruh"):
        super().__init__(ga, n, t, K, zk_general, zk_special, transform)
        self.sss = SubgroupShamir(ga.N, ga.factors, n, t, bound=2 * K + 1)
        self.mult = exceptional_multipliers(K, self.sss.M)

    def keygen(self):
        METER.set_phase("keygen")
        METER.set_party(0)
        ga = self.ga
        self.x = self.sss.secret()
        self.pk = [ga.act(ga.E0, self.sss.exponent(c * self.x)) for c in self.mult]
        sh = self.sss.share(self.x)
        self.shares = {i: sh[i] for i in range(1, self.n + 1)}
        METER.send(curves=self.K)

    def secret_bits_per_party(self):
        return self.sss.share_bits()

    def key_term(self, i, c, S):
        eff = self.sss.effective(i, self.shares[i], S)
        return (1 if c > 0 else -1) * self.mult[abs(c) - 1] * eff
