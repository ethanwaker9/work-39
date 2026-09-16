import math
import secrets
from itertools import combinations

from .galois import GaloisRing


def lagrange_at_zero(points, modulus):
    out = {}
    for i in points:
        num = 1
        den = 1
        for j in points:
            if j != i:
                num = num * j % modulus
                den = den * (j - i) % modulus
        out[i] = num * pow(den, -1, modulus) % modulus
    return out


class SubgroupShamir:
    name = "subgroup-shamir"

    def __init__(self, N, factors, n, t, bound=None):
        self.N = N
        self.n = n
        self.t = t
        limit = n if bound is None else max(n, bound)
        M = 1
        for p, e in factors.items():
            if p > limit:
                M *= p ** e
        self.M = M
        self.scale = N // M
        self.small = [p for p in factors if p <= limit]

    def secret(self):
        return secrets.randbelow(self.M)

    def exponent(self, x):
        return (self.scale * x) % self.N

    def share(self, s):
        coeffs = [s % self.M] + [secrets.randbelow(self.M) for _ in range(self.t - 1)]
        out = {}
        for i in range(1, self.n + 1):
            v = 0
            for c in reversed(coeffs):
                v = (v * i + c) % self.M
            out[i] = v
        return out

    def effective(self, i, share_i, S):
        L = lagrange_at_zero(sorted(S), self.M)
        return self.exponent(L[i] * share_i)

    def share_bits(self):
        return math.log2(self.M)

    def security_loss_bits(self):
        return math.log2(self.scale)


class LiftedShamir:
    name = "lifted-shamir"

    def __init__(self, N, factors, n, t):
        self.N = N
        self.n = n
        self.t = t
        self.components = []
        for p, e in sorted(factors.items()):
            q = p ** e
            d = 1
            while p ** d < n:
                d += 1
            R = GaloisRing(p, e, d)
            u = (N // q) * pow(N // q, -1, q) % N
            self.components.append((p, e, d, q, R, u))
        self._lag_cache = {}

    def secret(self):
        return secrets.randbelow(self.N)

    def exponent(self, x):
        return x % self.N

    def _is_inf(self, R, i):
        return i == R.p ** R.d

    def share(self, s):
        out = {i: [] for i in range(1, self.n + 1)}
        for (p, e, d, q, R, u) in self.components:
            coeffs = [R.scalar(s % q)] + [R.rand() for _ in range(self.t - 1)]
            for i in range(1, self.n + 1):
                if self._is_inf(R, i):
                    out[i].append(coeffs[-1])
                    continue
                x = R.point(i)
                v = R.zero()
                for c in reversed(coeffs):
                    v = R.add(R.mul(v, x), c)
                out[i].append(v)
        return out

    def lagrange(self, S):
        key = tuple(sorted(S))
        if key in self._lag_cache:
            return self._lag_cache[key]
        res = []
        for (p, e, d, q, R, u) in self.components:
            finite = [i for i in key if not self._is_inf(R, i)]
            pts = {i: R.point(i) for i in finite}
            lam = {}
            for i in finite:
                num = R.one()
                den = R.one()
                for j in finite:
                    if j != i:
                        num = R.mul(num, pts[j])
                        den = R.mul(den, R.sub(pts[j], pts[i]))
                lam[i] = R.mul(num, R.inv(den))
            for i in key:
                if self._is_inf(R, i):
                    prod = R.one()
                    for j in finite:
                        prod = R.mul(prod, R.sub(R.zero(), pts[j]))
                    lam[i] = prod
            res.append(lam)
        self._lag_cache[key] = res
        return res

    def effective(self, i, share_i, S):
        lams = self.lagrange(S)
        total = 0
        for idx, (p, e, d, q, R, u) in enumerate(self.components):
            v = R.mul(lams[idx][i], share_i[idx])
            total = (total + v[0] * u) % self.N
        return total

    def share_bits(self):
        return sum(R.bits() for (p, e, d, q, R, u) in self.components)

    def security_loss_bits(self):
        return 0.0


class Additive:
    name = "additive"

    def __init__(self, N, n):
        self.N = N
        self.n = n
        self.t = n

    def secret(self):
        return secrets.randbelow(self.N)

    def exponent(self, x):
        return x % self.N

    def share(self, s):
        parts = [secrets.randbelow(self.N) for _ in range(self.n - 1)]
        parts.append((s - sum(parts)) % self.N)
        return {i + 1: parts[i] for i in range(self.n)}

    def effective(self, i, share_i, S):
        return share_i % self.N

    def share_bits(self):
        return math.log2(self.N)


class Replicated:
    name = "replicated"

    def __init__(self, N, n, t):
        self.N = N
        self.n = n
        self.t = t
        self.maximal = list(combinations(range(1, n + 1), t - 1))

    def secret(self):
        return secrets.randbelow(self.N)

    def exponent(self, x):
        return x % self.N

    def share(self, s):
        parts = [secrets.randbelow(self.N) for _ in range(len(self.maximal) - 1)]
        parts.append((s - sum(parts)) % self.N)
        table = dict(zip(self.maximal, parts))
        out = {i: {} for i in range(1, self.n + 1)}
        for T, v in table.items():
            for i in range(1, self.n + 1):
                if i not in T:
                    out[i][T] = v
        return out

    def assign(self, S):
        S = sorted(S)
        return {T: min(i for i in S if i not in T) for T in self.maximal}

    def effective(self, i, share_i, S):
        owner = self.assign(S)
        return sum(v for T, v in share_i.items() if owner[T] == i) % self.N

    def shares_per_party(self):
        return math.comb(self.n - 1, self.t - 1)

    def share_bits(self):
        return self.shares_per_party() * math.log2(self.N)


class RecursiveSharing:
    name = "recursive"

    def __init__(self, N, n, t):
        self.N = N
        self.n = n
        self.t = t

    def secret(self):
        return secrets.randbelow(self.N)

    def exponent(self, x):
        return x % self.N

    def _keygen(self, g, parties, t, label, out):
        n = len(parties)
        if t == 1:
            for P in parties:
                out[P][label] = g
            return
        if t == n:
            parts = [secrets.randbelow(self.N) for _ in range(n - 1)]
            parts.append((g - sum(parts)) % self.N)
            for P, v in zip(parties, parts):
                out[P][label] = v
            return
        c = n // 2
        left = parties[:c]
        right = parties[c:]
        for ell in range(max(0, t - (n - c)), min(t, c) + 1):
            if ell == 0:
                self._keygen(g, right, t, label + (("R", t),), out)
            elif ell == t:
                self._keygen(g, left, t, label + (("L", t),), out)
            else:
                g1 = secrets.randbelow(self.N)
                g2 = (g - g1) % self.N
                self._keygen(g1, left, ell, label + (("L", ell, t),), out)
                self._keygen(g2, right, t - ell, label + (("R", t - ell, t),), out)

    def share(self, s):
        out = {i: {} for i in range(1, self.n + 1)}
        self._keygen(s, list(range(1, self.n + 1)), self.t, (), out)
        return out

    def _recover(self, P, T, parties, t, label):
        n = len(parties)
        if t == 1 or t == n:
            return label
        c = n // 2
        left = parties[:c]
        right = parties[c:]
        ell = len([x for x in T if x in left])
        if ell == 0:
            return self._recover(P, T, right, t, label + (("R", t),))
        if ell == t:
            return self._recover(P, T, left, t, label + (("L", t),))
        if P in left:
            return self._recover(P, [x for x in T if x in left], left, ell, label + (("L", ell, t),))
        return self._recover(P, [x for x in T if x in right], right, t - ell, label + (("R", t - ell, t),))

    def effective(self, i, share_i, S):
        S = sorted(S)
        label = self._recover(i, S, list(range(1, self.n + 1)), len(S), ())
        return share_i[label] % self.N

    def shares_per_party(self, shares=None):
        if shares is None:
            shares = self.share(0)
        return max(len(v) for v in shares.values())

    def share_bits(self, shares=None):
        return self.shares_per_party(shares) * math.log2(self.N)
