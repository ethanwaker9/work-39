import math
from fractions import Fraction


def reduce_form(a, b, c):
    while True:
        if c < a:
            a, b, c = c, -b, a
            continue
        if b > a or b <= -a:
            k = (a - b) // (2 * a)
            b_new = b + 2 * k * a
            c = (b_new * b_new - (b * b - 4 * a * c)) // (4 * a)
            b = b_new
            continue
        if a == c and b < 0:
            b = -b
            continue
        return a, b, c


def xgcd(a, b):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        return -a, -x0, -y0
    return a, x0, y0


def compose(f, g, D):
    a1, b1, c1 = f
    a2, b2, c2 = g
    s = (b1 + b2) // 2
    d1, u1, v1 = xgcd(a1, a2)
    d, x, y = xgcd(d1, s)
    v = x * v1
    w = y
    A = a1 * a2 // (d * d)
    B = b2 + 2 * (a2 // d) * (v * (s - b2) - w * c2)
    B %= 2 * A
    C = (B * B - D) // (4 * A)
    return reduce_form(A, B, C)


def form_pow(f, e, D, identity):
    if e < 0:
        f = (f[0], -f[1], f[2])
        f = reduce_form(*f)
        e = -e
    r = identity
    base = f
    while e:
        if e & 1:
            r = compose(r, base, D)
        base = compose(base, base, D)
        e >>= 1
    return r


def identity_form(D):
    if D % 4 == 0:
        return (1, 0, -D // 4)
    return (1, 1, (1 - D) // 4)


def prime_form(ell, p, D):
    return reduce_form(ell, 2, (p + 1) // ell)


def bsgs_order(f, D, bound):
    I = identity_form(D)
    m = math.isqrt(bound) + 1
    table = {}
    cur = I
    for j in range(m):
        if cur not in table:
            table[cur] = j
        cur = compose(cur, f, D)
    step = form_pow(f, -m, D, I)
    cur = step
    found = None
    for i in range(1, m + 2):
        if cur in table:
            found = i * m + table[cur]
            break
        cur = compose(cur, step, D)
    if found is None:
        raise ValueError("order not found")
    n = found
    for q in _prime_factors(n):
        while n % q == 0 and form_pow(f, n // q, D, I) == I:
            n //= q
    return n


def _prime_factors(n):
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


def factorize(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def bsgs_log(g, h, N, D):
    I = identity_form(D)
    m = math.isqrt(N) + 1
    table = {}
    cur = I
    for j in range(m):
        table.setdefault(cur, j)
        cur = compose(cur, g, D)
    step = form_pow(g, -m, D, I)
    cur = h
    for i in range(m + 1):
        if cur in table:
            return (i * m + table[cur]) % N
        cur = compose(cur, step, D)
    return None


def lll(basis, delta=Fraction(99, 100)):
    B = [list(v) for v in basis]
    n = len(B)

    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))

    def gram_schmidt():
        Bs = []
        mu = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = Fraction(dot(B[i], Bs[j])) / dot(Bs[j], Bs[j]) if dot(Bs[j], Bs[j]) else Fraction(0)
                v = [vi - mu[i][j] * bj for vi, bj in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu

    Bs, mu = gram_schmidt()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [x - q * y for x, y in zip(B[k], B[j])]
                Bs, mu = gram_schmidt()
        if dot(Bs[k], Bs[k]) >= (delta - mu[k][k - 1] ** 2) * dot(Bs[k - 1], Bs[k - 1]):
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu = gram_schmidt()
            k = max(k - 1, 1)
    return B


def babai(basis, target):
    n = len(basis)
    Bs = []
    for i in range(n):
        v = [Fraction(x) for x in basis[i]]
        for j in range(i):
            den = sum(x * x for x in Bs[j])
            mu = sum(Fraction(a) * b for a, b in zip(basis[i], Bs[j])) / den
            v = [vi - mu * bj for vi, bj in zip(v, Bs[j])]
        Bs.append(v)
    b = [Fraction(x) for x in target]
    coeffs = [0] * n
    for i in range(n - 1, -1, -1):
        den = sum(x * x for x in Bs[i])
        c = round(sum(x * y for x, y in zip(b, Bs[i])) / den)
        coeffs[i] = c
        b = [bi - c * vi for bi, vi in zip(b, basis[i])]
    close = [0] * len(target)
    for c, v in zip(coeffs, basis):
        if c:
            close = [x + c * y for x, y in zip(close, v)]
    return close
