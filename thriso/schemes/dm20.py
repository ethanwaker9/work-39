import secrets

from ..hashing import hash_to_int
from ..meter import METER
from ..sharing import SubgroupShamir, lagrange_at_zero
from .base import Scheme


class DM20(Scheme):
    name = "dm20"
    label = "DM20"
    active = False
    robust = False
    assumptions = "MT-GAIP, Power-DDHA"
    model = "ROM"

    def __init__(self, ga, n, t, K=16, **kw):
        super().__init__(ga, n, t, K)
        self.sss = SubgroupShamir(ga.N, ga.factors, n, t)

    def keygen(self):
        METER.set_phase("keygen")
        METER.set_party(0)
        ga = self.ga
        self.secrets = [self.sss.secret() for _ in range(self.K)]
        self.pk = [ga.act(ga.E0, self.sss.exponent(s)) for s in self.secrets]
        shares = [self.sss.share(s) for s in self.secrets]
        self.shares = {i: [shares[l][i] for l in range(self.K)] for i in range(1, self.n + 1)}
        METER.send(curves=self.K)

    def secret_bits_per_party(self):
        return self.K * self.sss.share_bits()

    def setup(self, S, cheat=None):
        return True

    def effective(self, i, l, S):
        return self.sss.effective(i, self.shares[i][l], S)

    def sign(self, msg, S, cheat=None):
        METER.set_phase("sign")
        ga = self.ga
        N = ga.N
        S = list(S)
        T = self.T
        X = [ga.E0] * T
        nonces = {}
        for i in S:
            METER.begin_round()
            METER.set_party(i)
            b = [ga.rand() for _ in range(T)]
            nonces[i] = b
            X = [ga.act(X[k], b[k]) for k in range(T)]
            METER.send(curves=T)
            METER.end_round()
        chal = self.fish.challenge(X, msg)
        z = [0] * T
        METER.begin_round()
        for i in S:
            METER.set_party(i)
            for k in range(T):
                c = chal[k]
                v = nonces[i][k]
                if c != 0:
                    v -= (1 if c > 0 else -1) * self.effective(i, abs(c) - 1, S)
                z[k] = (z[k] + v) % N
            METER.send(scalars=T)
        METER.end_round()
        return (chal, z)


class DEMS24(DM20):
    name = "dems24"
    label = "DEMS24"
    assumptions = "MT-GAIP, Power-DDHA"
    model = "QROM"

    def rerandomize(self, rho_seed):
        ga = self.ga
        M = self.sss.M
        self.rho = [hash_to_int("dems-rho", M, rho_seed, l) for l in range(self.K)]
        self.poly = [[self.rho[l]] + [hash_to_int("dems-g", M, rho_seed, l, j) for j in range(1, self.t)] for l in range(self.K)]
        METER.set_phase("sign")
        METER.set_party(0)
        self.pk_rand = [ga.act(self.pk[l], self.sss.exponent(self.rho[l])) for l in range(self.K)]

    def effective(self, i, l, S):
        M = self.sss.M
        F = 0
        for c in reversed(self.poly[l]):
            F = (F * i + c) % M
        L = lagrange_at_zero(sorted(S), M)
        return self.sss.exponent(L[i] * ((self.shares[i][l] + F) % M))

    def sign(self, msg, S, cheat=None):
        self.rerandomize(secrets.token_bytes(16))
        return super().sign(msg, S, cheat)

    def verify(self, msg, sig):
        METER.set_phase("verify")
        return self.fish.verify(self.pk_rand, msg, sig)
