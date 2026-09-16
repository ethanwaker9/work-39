import math

from ..hashing import challenge_vector, xof
from ..meter import METER

PARAMS = {
    1: {"K": 1, "T": 71, "slow": 16},
    16: {"K": 16, "T": 23, "slow": 15},
    256: {"K": 256, "T": 13, "slow": 12},
}


def exceptional_multipliers(K, N):
    out = [1]
    c = 2
    while len(out) < K:
        if all(math.gcd(c - x, N) == 1 and math.gcd(c + x, N) == 1 for x in out) and math.gcd(2 * c, N) == 1:
            out.append(c)
        c += 1
    return out


class CSIFiSh:
    def __init__(self, ga, K, T):
        self.ga = ga
        self.K = K
        self.T = T

    def base_curve(self, pk, c):
        if c == 0:
            return self.ga.E0
        if c > 0:
            return pk[c - 1]
        return self.ga.twist(pk[-c - 1])

    def challenge(self, curves, msg):
        return challenge_vector("csifish", self.T, self.K, list(curves), msg)

    def verify(self, pk, msg, sig):
        chal, z = sig
        if len(chal) != self.T or len(z) != self.T:
            return False
        curves = [self.ga.act(self.base_curve(pk, c), zi) for c, zi in zip(chal, z)]
        return self.challenge(curves, msg) == list(chal)

    def signature_bytes(self):
        return self.T * self.ga.scalar_bytes + math.ceil(self.T * math.log2(2 * self.K + 1) / 8)

    def pk_bytes(self):
        return self.K * self.ga.curve_bytes


class Scheme:
    name = "scheme"
    label = "scheme"
    active = False
    robust = False
    assumptions = ""
    model = "QROM"

    def __init__(self, ga, n, t, K=16, zk_general=112, zk_special=71, transform="unruh"):
        self.ga = ga
        self.n = n
        self.t = t
        self.K = K
        self.T = PARAMS[K]["T"]
        self.fish = CSIFiSh(ga, K, self.T)
        self.zk_general = zk_general
        self.zk_special = zk_special
        self.transform = transform

    def sign_signature_bytes(self):
        return self.fish.signature_bytes()

    def verify(self, msg, sig):
        METER.set_phase("verify")
        return self.fish.verify(self.pk, msg, sig)

    def secret_bits_per_party(self):
        raise NotImplementedError

    def pk_bytes(self):
        return self.fish.pk_bytes()
