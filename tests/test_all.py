import itertools
import os
import random
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "bench"))

from thriso.groupaction import MockAction, ToyCSIDHAction, CSIDH512_N
from thriso.meter import METER
from thriso.registry import SCHEMES, ORDER, parties_for, signer_set
from thriso.sharing import LiftedShamir, SubgroupShamir
from make_tables import sign_formula, setup_formula

FAC512 = {3: 1, 37: 1, 1407181: 1, 51593604295295867744293584889: 1,
          31599414504681995853008278745587832204909: 1}


class TestGroupAction(unittest.TestCase):
    def test_toy_action_is_compatible(self):
        ga = ToyCSIDHAction(count=False)
        rng = random.Random(7)
        a = rng.randrange(ga.N)
        b = rng.randrange(ga.N)
        left = ga.act(ga.act(ga.E0, a), b)
        right = ga.act(ga.E0, (a + b) % ga.N)
        self.assertEqual(left, right)
        self.assertEqual(ga.twist(ga.act(ga.E0, a)), ga.act(ga.E0, (-a) % ga.N))
        self.assertTrue(ga.cs.is_supersingular(left))


class TestLiftedSharing(unittest.TestCase):
    def check(self, N, factors, n, t, trials=3):
        sss = LiftedShamir(N, factors, n, t)
        rng = random.Random(n * 1000 + t)
        committees = list(itertools.combinations(range(1, n + 1), t))
        rng.shuffle(committees)
        for _ in range(trials):
            s = sss.secret()
            shares = sss.share(s)
            for S in committees[:6]:
                total = sum(sss.effective(i, shares[i], list(S)) for i in S) % N
                self.assertEqual(total, s % N)

    def test_many_parties_csidh512(self):
        self.check(CSIDH512_N, FAC512, 40, 3)
        self.check(CSIDH512_N, FAC512, 81, 2)

    def test_point_at_infinity(self):
        self.check(CSIDH512_N, FAC512, 9, 4)
        self.check(CSIDH512_N, FAC512, 3, 2)

    def test_toy_class_number(self):
        ga = ToyCSIDHAction(count=False)
        self.check(ga.N, ga.factors, 7, 4)

    def test_share_sizes(self):
        self.assertAlmostEqual(LiftedShamir(CSIDH512_N, FAC512, 100, 2).share_bits(), 268.69, places=2)
        self.assertAlmostEqual(SubgroupShamir(CSIDH512_N, FAC512, 100, 2).share_bits(), 250.34, places=2)

    def test_privacy_of_small_coalitions(self):
        N = 15
        factors = {3: 1, 5: 1}
        n, t = 5, 3
        sss = LiftedShamir(N, factors, n, t)
        self.assertEqual([c[2] for c in sss.components], [2, 1])
        for coalition in ((1, 5), (2, 4), (3, 5)):
            for idx, (p, e, d, q, R, u) in enumerate(sss.components):
                elements = list(itertools.product(range(q), repeat=d))
                for residue in range(q):
                    views = set()
                    for coeffs in itertools.product(elements, repeat=t - 1):
                        stream = iter(coeffs)
                        original = R.rand
                        R.rand = lambda: next(stream)
                        others = [comp for k, comp in enumerate(sss.components) if k != idx]
                        saved = [comp[4].rand for comp in others]
                        for comp in others:
                            comp[4].rand = comp[4].zero
                        try:
                            shares = sss.share(residue)
                        finally:
                            R.rand = original
                            for comp, fn in zip(others, saved):
                                comp[4].rand = fn
                        views.add(tuple(shares[i][idx] for i in coalition))
                    self.assertEqual(len(views), len(elements) ** (t - 1))


class TestSchemes(unittest.TestCase):
    def run_scheme(self, key, t, ga, K=16):
        n, t = parties_for(key, t)
        S = signer_set(key, n, t)
        METER.reset()
        sch = SCHEMES[key](ga, n=n, t=t, K=K)
        sch.keygen()
        sch.setup(S)
        sig = sch.sign(b"test message", S)
        self.assertIsNotNone(sig)
        self.assertTrue(sch.verify(b"test message", sig))
        self.assertFalse(sch.verify(b"other message", sig))
        summary = METER.summary(ga.curve_bytes, ga.scalar_bytes)
        tne = sum(1 for c in sig[0] if c != 0)
        expected = sign_formula(key, t, n, sch.T, K, 71, 112, tne)
        self.assertEqual(summary["sign"]["group_actions"], expected)
        self.assertEqual(summary.get("setup", {"group_actions": 0})["group_actions"], setup_formula(key, t, K, 112))

    def test_all_schemes_counting_backend(self):
        ga = MockAction()
        for key in ORDER:
            for t in (2, 3):
                with self.subTest(scheme=key, t=t):
                    self.run_scheme(key, t, ga)

    def test_real_backend_ours_and_dm20(self):
        ga = ToyCSIDHAction()
        for key in ("dm20", "ours"):
            with self.subTest(scheme=key):
                self.run_scheme(key, 2, ga)


class TestDeviations(unittest.TestCase):
    def test_audits_detect_deviations(self):
        ga = MockAction()
        t, n = 4, 7
        sch = SCHEMES["ours"](ga, n=n, t=t, K=16)
        sch.keygen()
        S = list(range(1, t + 1))
        self.assertIsNotNone(sch.setup(S))
        for pos in range(t):
            self.assertIsNone(sch.setup(S, {"party": S[pos], "type": "setup"}))
            sch.setup(S)
        for typ in ("chain", "response"):
            for pos in range(t):
                for rep in range(3):
                    sig = sch.sign(b"m", S, {"party": S[pos], "type": typ, "position": rep})
                    self.assertTrue(sig is None or not sch.verify(b"m", sig))


if __name__ == "__main__":
    unittest.main()
