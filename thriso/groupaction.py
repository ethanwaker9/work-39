import json
import math
import os
import secrets

from .csidh import CSIDH, CSIDH512_ELLS
from .classgroup import prime_form, bsgs_order, bsgs_log, lll, babai, factorize
from .meter import METER

TOY_ELLS = [3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 47]
CSIDH512_N = (3 * 37 * 1407181 * 51593604295295867744293584889
              * 31599414504681995853008278745587832204909)
CSIDH512_CURVE_BYTES = 64
CSIDH512_SCALAR_BYTES = (CSIDH512_N.bit_length() + 7) // 8

_HERE = os.path.dirname(os.path.abspath(__file__))
_CACHE = os.path.join(_HERE, "toy_params.json")


def _build_toy():
    cs = CSIDH(TOY_ELLS)
    p = cs.p
    D = -4 * p
    forms = [prime_form(l, p, D) for l in TOY_ELLS]
    N = bsgs_order(forms[0], D, int(math.isqrt(p) * math.log(p) * 2))
    logs = [bsgs_log(forms[0], f, N, D) for f in forms]
    n = len(TOY_ELLS)
    basis = [[N] + [0] * (n - 1)]
    for i in range(1, n):
        basis.append([(-logs[i]) % N] + [1 if j == i else 0 for j in range(1, n)])
    R = lll(basis)
    data = {"ells": TOY_ELLS, "p": p, "N": N, "logs": logs, "basis": R}
    with open(_CACHE, "w") as fh:
        json.dump(data, fh)
    return data


def toy_params():
    if os.path.exists(_CACHE):
        with open(_CACHE) as fh:
            return json.load(fh)
    return _build_toy()


class GroupAction:
    name = "abstract"

    def __init__(self, N, curve_bytes, scalar_bytes):
        self.N = N
        self.curve_bytes = curve_bytes
        self.scalar_bytes = scalar_bytes
        self.E0 = 0
        self.factors = factorize(N) if N < 1 << 64 else {3: 1, 37: 1, 1407181: 1,
                                                          51593604295295867744293584889: 1,
                                                          31599414504681995853008278745587832204909: 1}

    def rand(self):
        return secrets.randbelow(self.N)

    def act(self, E, a):
        raise NotImplementedError

    def twist(self, E):
        raise NotImplementedError

    def encode(self, E):
        return int(E).to_bytes(self.curve_bytes, "big")


class ToyCSIDHAction(GroupAction):
    name = "toy-csidh"

    def __init__(self, count=True):
        prm = toy_params()
        self.cs = CSIDH(prm["ells"])
        self.p = prm["p"]
        self.basis = prm["basis"]
        self.n = len(prm["ells"])
        self.count = count
        super().__init__(prm["N"], (self.p.bit_length() + 7) // 8, (prm["N"].bit_length() + 7) // 8)

    def vector(self, a):
        tgt = [a % self.N] + [0] * (self.n - 1)
        c = babai(self.basis, tgt)
        return [x - y for x, y in zip(tgt, c)]

    def act(self, E, a):
        if self.count:
            METER.count_ga()
        a %= self.N
        if a == 0:
            return E
        return self.cs.act(E, self.vector(a))

    def twist(self, E):
        return (-E) % self.p


class MockAction(GroupAction):
    name = "mock"

    def __init__(self, N=CSIDH512_N, curve_bytes=CSIDH512_CURVE_BYTES, scalar_bytes=CSIDH512_SCALAR_BYTES):
        super().__init__(N, curve_bytes, scalar_bytes)

    def act(self, E, a):
        METER.count_ga()
        return (E + a) % self.N

    def twist(self, E):
        return (-E) % self.N


class CSIDH512Timer:
    def __init__(self, bound=5):
        self.cs = CSIDH(CSIDH512_ELLS)
        self.bound = bound

    def one(self):
        e = [secrets.randbelow(2 * self.bound + 1) - self.bound for _ in CSIDH512_ELLS]
        return self.cs.act(0, e)
