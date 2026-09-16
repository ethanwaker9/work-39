import secrets

from ..hashing import xof, hash_to_int, commit, LAMBDA_BYTES
from ..meter import METER
from ..nizk import SigmaGAIP
from ..sharing import LiftedShamir
from .base import Scheme


class Albacore(Scheme):
    name = "albacore"
    label = "Ours"
    active = True
    robust = False
    assumptions = "MT-GAIP"
    model = "QROM"

    def __init__(self, ga, n, t, K=16, zk_general=112, zk_special=71, transform="unruh"):
        super().__init__(ga, n, t, K, zk_general, zk_special, transform)
        self.sss = LiftedShamir(ga.N, ga.factors, n, t)
        self.sigma = SigmaGAIP(ga, zk_general, zk_special, "unruh", force_general=True)
        self.committees = {}

    def keygen(self):
        METER.set_phase("keygen")
        METER.set_party(0)
        ga = self.ga
        self.secrets = [self.sss.secret() for _ in range(self.K)]
        self.pk = [ga.act(ga.E0, s) for s in self.secrets]
        shares = [self.sss.share(s) for s in self.secrets]
        self.shares = {i: [shares[l][i] for l in range(self.K)] for i in range(1, self.n + 1)}
        METER.send(curves=self.K)

    def secret_bits_per_party(self):
        return self.K * self.sss.share_bits()

    def committee_bits_per_party(self):
        return self.K * self.ga.N.bit_length()

    def setup(self, S, cheat=None):
        METER.set_phase("setup")
        ga = self.ga
        N = ga.N
        S = list(S)
        t = len(S)
        K = self.K
        mu = {i: [0] * K for i in S}
        METER.begin_round()
        for a in range(t):
            for b in range(a + 1, t):
                i, j = S[a], S[b]
                METER.set_party(i)
                for l in range(K):
                    r = secrets.randbelow(N)
                    mu[i][l] = (mu[i][l] + r) % N
                    mu[j][l] = (mu[j][l] - r) % N
                METER.send(scalars=K)
        METER.end_round()
        eff = {}
        for i in S:
            eff[i] = [(self.sss.effective(i, self.shares[i][l], S) + mu[i][l]) % N for l in range(K)]
        chains = [[ga.E0] for _ in range(K)]
        proofs = {}
        ctx = ("setup", tuple(S))
        for pos, i in enumerate(S):
            METER.begin_round()
            METER.set_party(i)
            for prev in S[:pos]:
                if prev == S[-1]:
                    continue
                for l in range(K):
                    pairs = [(chains[l][S.index(prev)], chains[l][S.index(prev) + 1])]
                    if not self.sigma.verify(pairs, proofs[(prev, l)], (ctx, prev, l)):
                        METER.end_round()
                        return None
            for l in range(K):
                w = eff[i][l]
                if cheat is not None and cheat.get("party") == i and cheat.get("type") == "setup":
                    w = (w + 1) % N
                nxt = ga.act(chains[l][-1], w)
                chains[l].append(nxt)
                METER.send(curves=1)
                if pos < t - 1:
                    proofs[(i, l)] = self.sigma.prove([(chains[l][-2], nxt)], w, (ctx, i, l))
            METER.end_round()
        for l in range(K):
            if chains[l][-1] != self.pk[l]:
                return None
        com = {"S": S, "eff": eff, "chains": chains}
        self.committees[tuple(S)] = com
        return com

    def sign(self, msg, S, cheat=None):
        METER.set_phase("sign")
        ga = self.ga
        N = ga.N
        S = list(S)
        t = len(S)
        T = self.T
        com = self.committees.get(tuple(S))
        if com is None:
            raise ValueError("committee setup missing")
        eff = com["eff"]
        chains = com["chains"]
        sid = secrets.token_bytes(16)
        salts = {}
        coms = {}
        X = [[ga.E0] for _ in range(T)]
        nonces = {}
        for pos, i in enumerate(S):
            METER.begin_round()
            METER.set_party(i)
            kappa = secrets.token_bytes(LAMBDA_BYTES)
            c, rnd = commit(("salt", sid, pos + 1, kappa))
            salts[i] = (kappa, rnd)
            coms[i] = c
            METER.send(hashbytes=LAMBDA_BYTES)
            b = [ga.rand() for _ in range(T)]
            nonces[i] = b
            for k in range(T):
                applied = b[k]
                if cheat is not None and cheat.get("party") == i and cheat.get("type") == "chain" and k == cheat.get("position", 0):
                    applied = (b[k] + 1) % N
                X[k].append(ga.act(X[k][-1], applied))
            METER.send(curves=T)
            METER.end_round()
        METER.begin_round()
        for pos, i in enumerate(S):
            METER.set_party(i)
            kappa, rnd = salts[i]
            if commit(("salt", sid, pos + 1, kappa), rnd)[0] != coms[i]:
                METER.end_round()
                return None
            METER.send(hashbytes=len(kappa) + len(rnd))
        transcript = xof("albacore-transcript", sid, msg, S, [row for row in X], [salts[i][0] for i in S])
        METER.set_party(S[-1])
        rho = [hash_to_int("albacore-rho", N, transcript, k) for k in range(T)]
        Fstar = [ga.act(X[k][-1], rho[k]) for k in range(T)]
        METER.end_round()
        chal = self.fish.challenge(Fstar, msg)
        u = [0] * T
        for pos, i in enumerate(S):
            METER.begin_round()
            METER.set_party(i)
            if pos > 0:
                for k in range(T):
                    c = chal[k]
                    if c == 0:
                        continue
                    base = chains[abs(c) - 1][pos]
                    if c < 0:
                        base = ga.twist(base)
                    if ga.act(base, u[k]) != X[k][pos]:
                        METER.end_round()
                        return None
            for k in range(T):
                c = chal[k]
                v = nonces[i][k]
                if c != 0:
                    sgn = 1 if c > 0 else -1
                    v = v - sgn * eff[i][abs(c) - 1]
                if cheat is not None and cheat.get("party") == i and cheat.get("type") == "response" and k == cheat.get("position", 0):
                    v += 1
                u[k] = (u[k] + v) % N
            METER.send(scalars=T)
            METER.end_round()
        z = [(u[k] + rho[k]) % N for k in range(T)]
        return (chal, z)
