import math
import secrets


class CSIDH:
    def __init__(self, ells, p=None):
        self.ells = list(ells)
        self.p = p if p is not None else 4 * math.prod(self.ells) - 1
        self.n = len(self.ells)
        self.cofactor = (self.p + 1) // math.prod(self.ells)

    def is_square(self, v):
        v %= self.p
        if v == 0:
            return 0
        return 1 if pow(v, (self.p - 1) // 2, self.p) == 1 else -1

    def xdbl(self, X, Z, A24, C24):
        p = self.p
        t0 = (X - Z) % p
        t1 = (X + Z) % p
        t0 = t0 * t0 % p
        t1 = t1 * t1 % p
        Z2 = C24 * t0 % p
        X2 = Z2 * t1 % p
        t1 = (t1 - t0) % p
        t0 = A24 * t1 % p
        Z2 = (Z2 + t0) * t1 % p
        return X2, Z2

    def xadd(self, XP, ZP, XQ, ZQ, XD, ZD):
        p = self.p
        t0 = (XP + ZP) % p
        t1 = (XP - ZP) % p
        t2 = (XQ - ZQ) % p
        t3 = (XQ + ZQ) % p
        t0 = t0 * t2 % p
        t1 = t1 * t3 % p
        t2 = (t0 + t1) % p
        t3 = (t0 - t1) % p
        t2 = t2 * t2 % p
        t3 = t3 * t3 % p
        return ZD * t2 % p, XD * t3 % p

    def xmul(self, X, Z, k, A24, C24):
        if k == 0:
            return 1, 0
        if k == 1:
            return X, Z
        X0, Z0 = X, Z
        X1, Z1 = self.xdbl(X, Z, A24, C24)
        for bit in bin(k)[3:]:
            if bit == '1':
                X0, Z0 = self.xadd(X1, Z1, X0, Z0, X, Z)
                X1, Z1 = self.xdbl(X1, Z1, A24, C24)
            else:
                X1, Z1 = self.xadd(X0, Z0, X1, Z1, X, Z)
                X0, Z0 = self.xdbl(X0, Z0, A24, C24)
        return X0, Z0

    def isogeny(self, A, C, XK, ZK, ell, pushes):
        p = self.p
        A24 = (A + 2 * C) % p
        C24 = 4 * C % p
        k = (ell - 1) // 2
        a = A24
        d = (A - 2 * C) % p
        pp = 1
        pm = 1
        Xi, Zi = XK, ZK
        mults = []
        prevX, prevZ = None, None
        for i in range(1, k + 1):
            if i == 1:
                Xi, Zi = XK, ZK
            elif i == 2:
                prevX, prevZ = XK, ZK
                Xi, Zi = self.xdbl(XK, ZK, A24, C24)
            else:
                nX, nZ = self.xadd(Xi, Zi, XK, ZK, prevX, prevZ)
                prevX, prevZ = Xi, Zi
                Xi, Zi = nX, nZ
            mults.append((Xi, Zi))
            pp = pp * ((Xi + Zi) % p) % p
            pm = pm * ((Xi - Zi) % p) % p
        pp2 = pp * pp % p
        pp4 = pp2 * pp2 % p
        pp8 = pp4 * pp4 % p
        pm2 = pm * pm % p
        pm4 = pm2 * pm2 % p
        pm8 = pm4 * pm4 % p
        an = pow(a, ell, p) * pp8 % p
        dn = pow(d, ell, p) * pm8 % p
        An = 2 * (an + dn) % p
        Cn = (an - dn) % p
        out = []
        for (XQ, ZQ) in pushes:
            num = 1
            den = 1
            for (Xm, Zm) in mults:
                t0 = (XQ * Xm - ZQ * Zm) % p
                t1 = (XQ * Zm - ZQ * Xm) % p
                num = num * t0 % p
                den = den * t1 % p
            out.append((XQ * num % p * num % p, ZQ * den % p * den % p))
        return An, Cn, out

    def normalize(self, A, C):
        return A * pow(C, -1, self.p) % self.p

    def act(self, A, exps):
        p = self.p
        e = list(exps)
        A = A % p
        C = 1
        while any(e):
            x = secrets.randbelow(p - 2) + 2
            Ac = A * pow(C, -1, p) % p
            s = self.is_square(x * x % p * x % p + Ac * x % p * x % p + x)
            if s == 0:
                continue
            S = [i for i in range(self.n) if e[i] != 0 and (1 if e[i] > 0 else -1) == s]
            if not S:
                continue
            k = 1
            for i in S:
                k *= self.ells[i]
            A24 = (A + 2 * C) % p
            C24 = 4 * C % p
            XQ, ZQ = self.xmul(x, 1, (p + 1) // k, A24, C24)
            for idx in reversed(S):
                ell = self.ells[idx]
                A24 = (A + 2 * C) % p
                C24 = 4 * C % p
                XR, ZR = self.xmul(XQ, ZQ, k // ell, A24, C24)
                if ZR != 0:
                    A, C, pushed = self.isogeny(A, C, XR, ZR, ell, [(XQ, ZQ)])
                    XQ, ZQ = pushed[0]
                    e[idx] -= s
                k //= ell
            A = A * pow(C, -1, p) % p
            C = 1
        return A % p

    def is_supersingular(self, A, trials=2):
        p = self.p
        for _ in range(trials):
            x = secrets.randbelow(p - 2) + 2
            A24 = (A + 2) % p
            C24 = 4
            XQ, ZQ = self.xmul(x, 1, 4, A24, C24)
            order_part = 1
            for ell in self.ells:
                XR, ZR = self.xmul(XQ, ZQ, (p + 1) // (4 * ell), A24, C24)
                if ZR != 0:
                    XT, ZT = self.xmul(XR, ZR, ell, A24, C24)
                    if ZT != 0:
                        return False
                    order_part *= ell
                    if order_part > 4 * math.isqrt(p):
                        return True
        return False

    def twist(self, A):
        return (-A) % self.p


def first_odd_primes(k, skip=()):
    out = []
    c = 3
    while len(out) < k:
        if all(c % q for q in range(3, int(c ** 0.5) + 1, 2)) and c not in skip:
            out.append(c)
        c += 2
    return out


CSIDH512_ELLS = first_odd_primes(73) + [587]
