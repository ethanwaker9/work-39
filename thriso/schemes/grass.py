import math
import secrets
from itertools import combinations

from ..hashing import commit, challenge_vector, xof, LAMBDA_BYTES
from ..meter import METER
from ..nizk import SigmaGAIP
from .base import Scheme


def salted_challenge(fish, curves, salt, msg):
    return challenge_vector("csifish-salted", fish.T, fish.K, list(curves), salt, msg)


class GRASS(Scheme):
    name = "grass"
    label = "BBMP24"
    active = True
    robust = False
    assumptions = "MT-GAIP, 2-weak pseudorandomness"
    model = "QROM"

    def __init__(self, ga, n, t, K=16, zk_general=112, zk_special=71, transform="unruh"):
        super().__init__(ga, n, t, K, zk_general, zk_special, transform)
        self.sigma = SigmaGAIP(ga, zk_general, zk_special, transform)
        maximal = list(combinations(range(1, n + 1), t - 1))
        self.shards = [tuple(i for i in range(1, n + 1) if i not in U) for U in maximal]
        self.M = len(self.shards)

    def keygen(self):
        METER.set_phase("keygen")
        ga = self.ga
        N = ga.N
        self.shard_values = []
        self.inter = []
        self.pk = []
        for l in range(self.K):
            vals = []
            chain = [ga.E0]
            for m, shard in enumerate(self.shards):
                owner = shard[0]
                METER.begin_round()
                METER.set_party(owner)
                v = ga.rand()
                vals.append(v)
                nxt = ga.act(chain[-1], v)
                proof = self.sigma.prove([(chain[-1], nxt)], v, ("grass-kg", l, m))
                METER.send(curves=1)
                METER.end_round()
                METER.begin_round()
                for j in range(1, self.n + 1):
                    if j == owner:
                        continue
                    METER.set_party(j)
                    if not self.sigma.verify([(chain[-1], nxt)], proof, ("grass-kg", l, m)):
                        raise ValueError("keygen proof")
                METER.end_round()
                chain.append(nxt)
            self.shard_values.append(vals)
            self.inter.append(chain)
            self.pk.append(chain[-1])

    def secret_bits_per_party(self):
        per = sum(1 for sh in self.shards if 1 in sh)
        return self.K * per * math.log2(self.ga.N)

    def setup(self, S, cheat=None):
        return True

    def inter_curve(self, c, m):
        if c == 0:
            return self.ga.E0
        base = self.inter[abs(c) - 1][m]
        return base if c > 0 else self.ga.twist(base)

    def sign(self, msg, S, cheat=None):
        METER.set_phase("sign")
        ga = self.ga
        N = ga.N
        S = sorted(S)
        T = self.T
        tau = [min(i for i in shard if i in S) for shard in self.shards]
        METER.begin_round()
        salts = {}
        coms = {}
        for i in S:
            METER.set_party(i)
            salts[i] = secrets.token_bytes(LAMBDA_BYTES)
            coms[i] = commit(("salt", i, salts[i]))
            METER.send(hashbytes=LAMBDA_BYTES)
        METER.end_round()
        X = [[ga.E0] for _ in range(T)]
        nonces = []
        for m in range(self.M):
            METER.begin_round()
            METER.set_party(tau[m])
            g = [ga.rand() for _ in range(T)]
            nonces.append(g)
            for k in range(T):
                applied = g[k]
                if cheat is not None and cheat.get("party") == tau[m] and cheat.get("type") == "chain" and k == cheat.get("position", 0) and m == S.index(tau[m]):
                    applied = (g[k] + 1) % N
                X[k].append(ga.act(X[k][-1], applied))
            METER.send(curves=T)
            METER.end_round()
        METER.begin_round()
        for i in S:
            METER.set_party(i)
            if commit(("salt", i, salts[i]), coms[i][1])[0] != coms[i][0]:
                METER.end_round()
                return None
            METER.send(hashbytes=2 * LAMBDA_BYTES)
        METER.end_round()
        salt = xof("grass-salt", [salts[i] for i in S])
        chal = salted_challenge(self.fish, [X[k][-1] for k in range(T)], salt, msg)
        u = [0] * T
        for m in range(self.M):
            METER.begin_round()
            METER.set_party(tau[m])
            for k in range(T):
                c = chal[k]
                v = nonces[m][k]
                if c != 0:
                    v -= (1 if c > 0 else -1) * self.shard_values[abs(c) - 1][m]
                u[k] = (u[k] + v) % N
            METER.send(scalars=T)
            for j in S:
                METER.set_party(j)
                for k in range(T):
                    if ga.act(self.inter_curve(chal[k], m + 1), u[k]) != X[k][m + 1]:
                        METER.end_round()
                        return None
            METER.end_round()
        return (chal, u, salt)

    def verify(self, msg, sig):
        METER.set_phase("verify")
        chal, z, salt = sig
        curves = [self.ga.act(self.fish.base_curve(self.pk, c), zi) for c, zi in zip(chal, z)]
        return salted_challenge(self.fish, curves, salt, msg) == list(chal)

    def sign_signature_bytes(self):
        return self.fish.signature_bytes() + LAMBDA_BYTES


class GRASSPlus(GRASS):
    name = "grassplus"
    label = "BBDMP25"
    assumptions = "Chain-GAIP, Graph-GAIP"
    model = "ROM, BBGAM"

    def __init__(self, ga, n, t, K=16, zk_general=112, zk_special=71, transform="unruh"):
        Scheme.__init__(self, ga, n, t, K, zk_general, zk_special, transform)
        self.sigma = SigmaGAIP(ga, zk_general, zk_special, "unruh")

    def _rec(self, l, g, x, parties, t, label, out):
        ga = self.ga
        N = ga.N
        n = len(parties)
        if t == 1:
            y = ga.act(x, g)
            for P in parties:
                out[P][label] = (g, x, y)
            return
        if t == n:
            parts = [secrets.randbelow(N) for _ in range(n - 1)]
            parts.append((g - sum(parts)) % N)
            cur = x
            for P, v in zip(parties, parts):
                nxt = ga.act(cur, v)
                out[P][label] = (v, cur, nxt)
                cur = nxt
            return
        c = n // 2
        left = parties[:c]
        right = parties[c:]
        for ell in range(max(0, t - (n - c)), min(t, c) + 1):
            if ell == 0:
                self._rec(l, g, x, right, t, label + (("R", t),), out)
            elif ell == t:
                self._rec(l, g, x, left, t, label + (("L", t),), out)
            else:
                g1 = secrets.randbelow(N)
                x1 = ga.act(x, g1)
                self._rec(l, g1, x, left, ell, label + (("L", ell, t),), out)
                self._rec(l, (g - g1) % N, x1, right, t - ell, label + (("R", t - ell, t),), out)

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

    def keygen(self):
        METER.set_phase("keygen")
        METER.set_party(0)
        ga = self.ga
        self.store = []
        self.pk = []
        for l in range(self.K):
            s = ga.rand()
            out = {i: {} for i in range(1, self.n + 1)}
            self._rec(l, s, ga.E0, list(range(1, self.n + 1)), self.t, (), out)
            self.store.append(out)
            self.pk.append(ga.act(ga.E0, s))

    def secret_bits_per_party(self):
        per = max(len(self.store[0][i]) for i in self.store[0])
        return self.K * per * math.log2(self.ga.N)

    def sign(self, msg, S, cheat=None):
        METER.set_phase("sign")
        ga = self.ga
        N = ga.N
        S = sorted(S)
        T = self.T
        t = len(S)
        keys = {i: [self.store[l][i][self._recover(i, S, list(range(1, self.n + 1)), t, ())] for l in range(self.K)] for i in S}
        X = [[ga.E0] for _ in range(T)]
        nonces = {}
        proofs = {}
        salts = {}
        coms = {}
        for pos, i in enumerate(S):
            METER.begin_round()
            METER.set_party(i)
            for j in S[:pos]:
                for k in range(T):
                    if not self.sigma.verify([(X[k][S.index(j)], X[k][S.index(j) + 1])], proofs[j][k], ("gp", j, k)):
                        METER.end_round()
                        return None
            salts[i] = secrets.token_bytes(LAMBDA_BYTES)
            coms[i] = commit(("salt", i, salts[i]))
            g = [ga.rand() for _ in range(T)]
            nonces[i] = g
            pr = []
            for k in range(T):
                applied = g[k]
                if cheat is not None and cheat.get("party") == i and cheat.get("type") == "chain" and k == cheat.get("position", 0):
                    applied = (g[k] + 1) % N
                nxt = ga.act(X[k][-1], applied)
                pr.append(self.sigma.prove([(X[k][-1], nxt)], g[k], ("gp", i, k)))
                X[k].append(nxt)
            proofs[i] = pr
            METER.send(curves=T, hashbytes=LAMBDA_BYTES)
            METER.end_round()
        METER.begin_round()
        for i in S:
            METER.set_party(i)
            METER.send(hashbytes=2 * LAMBDA_BYTES)
        METER.end_round()
        salt = xof("grass-salt", [salts[i] for i in S])
        chal = salted_challenge(self.fish, [X[k][-1] for k in range(T)], salt, msg)
        rsp = {}
        for pos, i in enumerate(S):
            METER.begin_round()
            METER.set_party(i)
            for j in S:
                if j == i:
                    continue
                if commit(("salt", j, salts[j]), coms[j][1])[0] != coms[j][0]:
                    METER.end_round()
                    return None
                for k in range(T):
                    if not self.sigma.verify([(X[k][S.index(j)], X[k][S.index(j) + 1])], proofs[j][k], ("gp", j, k)):
                        METER.end_round()
                        return None
            for j in S[:pos]:
                jp = S.index(j)
                for k in range(T):
                    c = chal[k]
                    if c == 0:
                        ok = ga.act(X[k][jp], rsp[j][k]) == X[k][jp + 1]
                    else:
                        base = keys[j][abs(c) - 1][2]
                        base = base if c > 0 else ga.twist(base)
                        ok = ga.act(base, rsp[j][k]) == X[k][jp + 1]
                    if not ok:
                        METER.end_round()
                        return None
            out = []
            for k in range(T):
                c = chal[k]
                if c == 0:
                    out.append(nonces[i][k])
                else:
                    prev = rsp[S[pos - 1]][k] if pos > 0 else 0
                    out.append((prev + nonces[i][k] - (1 if c > 0 else -1) * keys[i][abs(c) - 1][0]) % N)
            rsp[i] = out
            METER.send(scalars=T)
            METER.end_round()
        z = []
        for k in range(T):
            if chal[k] == 0:
                z.append(sum(rsp[i][k] for i in S) % N)
            else:
                z.append(rsp[S[-1]][k])
        return (chal, z, salt)
