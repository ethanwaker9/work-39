import math
import secrets

from .hashing import xof, challenge_vector, commit, LAMBDA_BYTES
from .meter import METER


class SigmaGAIP:
    def __init__(self, ga, rounds_general=112, rounds_special=71, transform="unruh", force_general=False):
        self.ga = ga
        self.rounds_general = rounds_general
        self.rounds_special = rounds_special
        self.transform = transform
        self.force_general = force_general

    def _special(self, pairs):
        if self.force_general:
            return False
        return all(X == self.ga.E0 for X, _ in pairs)

    def _base(self, pairs, j, c):
        X, Y = pairs[j]
        if c == 0:
            return X
        if c == 1:
            return Y
        return self.ga.twist(Y)

    def rounds(self, pairs):
        return self.rounds_special if self._special(pairs) else self.rounds_general

    def prove(self, pairs, w, ctx):
        ga = self.ga
        special = self._special(pairs)
        K = 1
        R = self.rounds_special if special else self.rounds_general
        cset = [-1, 0, 1] if special else [0, 1]
        bs = [ga.rand() for _ in range(R)]
        coms = [[ga.act(X, b) for (X, _) in pairs] for b in bs]
        if self.transform == "fs":
            if special:
                chal = challenge_vector("sigma-fs", R, K, ctx, pairs, coms)
            else:
                chal = [c & 1 for c in _bits(xof("sigma-fs", ctx, pairs, coms, length=(R + 7) // 8), R)]
            resp = [(b - c * w) % ga.N for b, c in zip(bs, chal)]
            METER.send(scalars=R, hashbytes=(R * (2 if special else 1) + 7) // 8)
            return ("fs", chal, resp)
        hashes = []
        allresp = []
        for r, b in enumerate(bs):
            row = {}
            hrow = {}
            for c in cset:
                v = (b - c * w) % ga.N
                row[c] = v
                hrow[c] = xof("unruh-g", ctx, r, c, v, length=LAMBDA_BYTES)
            allresp.append(row)
            hashes.append(hrow)
        hlist = [[hashes[r][c] for c in cset] for r in range(R)]
        if special:
            chal = challenge_vector("sigma-unruh", R, K, ctx, pairs, coms, hlist)
        else:
            chal = [c & 1 for c in _bits(xof("sigma-unruh", ctx, pairs, coms, hlist, length=(R + 7) // 8), R)]
        opened = [allresp[r][chal[r]] for r in range(R)]
        unopened = [[hashes[r][c] for c in cset if c != chal[r]] for r in range(R)]
        METER.send(scalars=R, hashbytes=R * (len(cset) - 1) * LAMBDA_BYTES + (R * (2 if special else 1) + 7) // 8)
        return ("unruh", chal, opened, unopened)

    def verify(self, pairs, proof, ctx):
        ga = self.ga
        special = self._special(pairs)
        K = 1
        R = self.rounds_special if special else self.rounds_general
        cset = [-1, 0, 1] if special else [0, 1]
        kind = proof[0]
        chal = proof[1]
        if len(chal) != R:
            return False
        if kind == "fs":
            resp = proof[2]
            coms = [[ga.act(self._base(pairs, j, c), v) for j in range(len(pairs))] for v, c in zip(resp, chal)]
            if special:
                chal2 = challenge_vector("sigma-fs", R, K, ctx, pairs, coms)
            else:
                chal2 = [c & 1 for c in _bits(xof("sigma-fs", ctx, pairs, coms, length=(R + 7) // 8), R)]
            return chal2 == list(chal)
        opened = proof[2]
        unopened = proof[3]
        coms = [[ga.act(self._base(pairs, j, c), v) for j in range(len(pairs))] for v, c in zip(opened, chal)]
        hlist = []
        for r in range(R):
            others = list(unopened[r])
            row = []
            for c in cset:
                if c == chal[r]:
                    row.append(xof("unruh-g", ctx, r, c, opened[r], length=LAMBDA_BYTES))
                else:
                    row.append(others.pop(0))
            hlist.append(row)
        if special:
            chal2 = challenge_vector("sigma-unruh", R, K, ctx, pairs, coms, hlist)
        else:
            chal2 = [c & 1 for c in _bits(xof("sigma-unruh", ctx, pairs, coms, hlist, length=(R + 7) // 8), R)]
        return chal2 == list(chal)


def _bits(raw, R):
    v = int.from_bytes(raw, "big")
    return [(v >> i) & 1 for i in range(R)]


def poly_eval(coeffs, x, N):
    v = 0
    for c in reversed(coeffs):
        v = (v * x + c) % N
    return v


class PVP:
    def __init__(self, ga, reps=112, modulus=None):
        self.ga = ga
        self.reps = reps
        self.modulus = modulus if modulus is not None else ga.N

    def prove(self, f, base_pairs, xs, ctx, scale=1):
        ga = self.ga
        N = self.modulus
        t = len(f) - 1
        bs = [[secrets.randbelow(N) for _ in range(t + 1)] for _ in range(self.reps)]
        Fhat = [[ga.act(F, (scale * b[0]) % ga.N) for (F, _) in base_pairs] for b in bs]
        y0 = secrets.token_bytes(LAMBDA_BYTES)
        y0p = secrets.token_bytes(LAMBDA_BYTES)
        C0, _ = commit(("fhat", Fhat), y0)
        C0p, _ = commit(("stmt", base_pairs), y0p)
        C = {0: C0}
        Cp = {0: C0p}
        pieces = {0: (y0, y0p)}
        for i in xs:
            yi = secrets.token_bytes(LAMBDA_BYTES)
            yip = secrets.token_bytes(LAMBDA_BYTES)
            C[i], _ = commit(("b", [poly_eval(b, i, N) for b in bs]), yi)
            Cp[i], _ = commit(("x", xs[i]), yip)
            pieces[i] = (yi, yip)
        keys = sorted(C)
        d = _bits(xof("pvp", ctx, [C[k] for k in keys], [Cp[k] for k in keys], length=(self.reps + 7) // 8), self.reps)
        r = [[(bj - dj * fj) % N for bj, fj in zip(b, f)] for b, dj in zip(bs, d)]
        main = (C, Cp, r, d)
        METER.send(scalars=self.reps * (t + 1), hashbytes=2 * LAMBDA_BYTES * len(keys) + 2 * LAMBDA_BYTES * len(keys))
        return main, pieces

    def verify_piece(self, i, x, main, piece, base_pairs, ctx, scale=1):
        ga = self.ga
        N = self.modulus
        C, Cp, r, d = main
        keys = sorted(C)
        d2 = _bits(xof("pvp", ctx, [C[k] for k in keys], [Cp[k] for k in keys], length=(self.reps + 7) // 8), self.reps)
        if d2 != list(d):
            return False
        y, yp = piece
        if i == 0:
            if commit(("stmt", base_pairs), yp)[0] != Cp[0]:
                return False
            Fhat = []
            for rj, dj in zip(r, d):
                row = []
                for (F, Fp) in base_pairs:
                    row.append(ga.act(Fp if dj else F, (scale * rj[0]) % ga.N))
                Fhat.append(row)
            return commit(("fhat", Fhat), y)[0] == C[0]
        if commit(("x", x), yp)[0] != Cp[i]:
            return False
        vals = [(poly_eval(rj, i, N) + dj * x) % N for rj, dj in zip(r, d)]
        return commit(("b", vals), y)[0] == C[i]


class NITZKShamir:
    def __init__(self, N, lam=128, exceptional=None):
        self.N = N
        self.exceptional = exceptional if exceptional is not None else N
        self.reps = max(1, math.ceil(lam / math.log2(self.exceptional)))

    def prove(self, f, xs, ctx):
        N = self.N
        t = len(f) - 1
        bs = [[secrets.randbelow(N) for _ in range(t + 1)] for _ in range(self.reps)]
        C = {}
        Cp = {}
        pieces = {}
        for i in xs:
            yi = secrets.token_bytes(LAMBDA_BYTES)
            yip = secrets.token_bytes(LAMBDA_BYTES)
            C[i], _ = commit(("b", [poly_eval(b, i, N) for b in bs]), yi)
            Cp[i], _ = commit(("x", xs[i]), yip)
            pieces[i] = (yi, yip)
        keys = sorted(C)
        dvec = [xof("nitzk", ctx, k, [C[q] for q in keys], [Cp[q] for q in keys], length=8) for k in range(self.reps)]
        ds = [int.from_bytes(v, "big") % self.exceptional for v in dvec]
        r = [[(bj - dj * fj) % N for bj, fj in zip(b, f)] for b, dj in zip(bs, ds)]
        METER.send(scalars=self.reps * (t + 1), hashbytes=4 * LAMBDA_BYTES * len(keys))
        return (C, Cp, r), pieces

    def verify(self, i, x, main, piece, ctx):
        N = self.N
        C, Cp, r = main
        keys = sorted(C)
        dvec = [xof("nitzk", ctx, k, [C[q] for q in keys], [Cp[q] for q in keys], length=8) for k in range(self.reps)]
        ds = [int.from_bytes(v, "big") % self.exceptional for v in dvec]
        y, yp = piece
        if commit(("x", x), yp)[0] != Cp[i]:
            return False
        vals = [(poly_eval(rj, i, N) + dj * x) % N for rj, dj in zip(r, ds)]
        return commit(("b", vals), y)[0] == C[i]
