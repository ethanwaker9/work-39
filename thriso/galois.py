import secrets


def _poly_mod_p(a, f, p):
    a = [x % p for x in a]
    df = len(f) - 1
    while len(a) - 1 >= df and any(a):
        while a and a[-1] == 0:
            a.pop()
        if len(a) - 1 < df:
            break
        c = a[-1]
        shift = len(a) - 1 - df
        for i in range(df + 1):
            a[shift + i] = (a[shift + i] - c * f[i]) % p
        while a and a[-1] == 0:
            a.pop()
    return a


def _poly_mul_p(a, b, p):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
    return out


def _poly_gcd_p(a, b, p):
    a = [x % p for x in a]
    b = [x % p for x in b]
    while a and a[-1] == 0:
        a.pop()
    while b and b[-1] == 0:
        b.pop()
    while b:
        inv = pow(b[-1], -1, p)
        bm = [(x * inv) % p for x in b]
        a = _poly_mod_p(a, bm, p)
        a, b = b, a
    return a


def _xpow_mod(e, f, p):
    result = [1]
    base = [0, 1]
    while e:
        if e & 1:
            result = _poly_mod_p(_poly_mul_p(result, base, p), f, p)
        base = _poly_mod_p(_poly_mul_p(base, base, p), f, p)
        e >>= 1
    return result


def _prime_divisors(n):
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def _sub_x(h, p):
    h = list(h) + [0] * max(0, 2 - len(h))
    h[1] = (h[1] - 1) % p
    while h and h[-1] == 0:
        h.pop()
    return h


def is_irreducible(f, p):
    d = len(f) - 1
    if d == 1:
        return True
    if _sub_x(_xpow_mod(p ** d, f, p), p):
        return False
    for q in _prime_divisors(d):
        g = _poly_gcd_p(f, _sub_x(_xpow_mod(p ** (d // q), f, p), p), p)
        if len(g) > 1:
            return False
    return True


def find_irreducible(p, d):
    if d == 1:
        return [0, 1]
    c = 0
    while True:
        coeffs = []
        x = c
        for _ in range(d):
            coeffs.append(x % p)
            x //= p
        f = coeffs + [1]
        if f[0] != 0 and is_irreducible(f, p):
            return f
        c += 1


class GaloisRing:
    def __init__(self, p, e, d):
        self.p = p
        self.e = e
        self.d = d
        self.q = p ** e
        self.f = find_irreducible(p, d)
        self.unit_order = (p ** d - 1) * p ** ((e - 1) * d)

    def zero(self):
        return (0,) * self.d

    def one(self):
        return (1,) + (0,) * (self.d - 1)

    def scalar(self, s):
        return (s % self.q,) + (0,) * (self.d - 1)

    def add(self, a, b):
        return tuple((x + y) % self.q for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple((x - y) % self.q for x, y in zip(a, b))

    def mul(self, a, b):
        d = self.d
        q = self.q
        if d == 1:
            return ((a[0] * b[0]) % q,)
        prod = [0] * (2 * d - 1)
        for i in range(d):
            ai = a[i]
            if ai:
                for j in range(d):
                    prod[i + j] += ai * b[j]
        f = self.f
        for k in range(2 * d - 2, d - 1, -1):
            c = prod[k] % q
            if c:
                for i in range(d):
                    prod[k - d + i] -= c * f[i]
            prod[k] = 0
        return tuple(x % q for x in prod[:d])

    def power(self, a, n):
        r = self.one()
        b = a
        while n:
            if n & 1:
                r = self.mul(r, b)
            b = self.mul(b, b)
            n >>= 1
        return r

    def inv(self, a):
        return self.power(a, self.unit_order - 1)

    def rand(self):
        return tuple(secrets.randbelow(self.q) for _ in range(self.d))

    def point(self, i):
        out = []
        x = i
        for _ in range(self.d):
            out.append(x % self.p)
            x //= self.p
        if x != 0:
            raise ValueError("exceptional set too small")
        return tuple(out)

    def bits(self):
        import math
        return self.d * self.e * math.log2(self.p)
