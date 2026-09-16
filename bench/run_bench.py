import argparse
import json
import math
import os
import random
import statistics
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thriso.groupaction import MockAction, ToyCSIDHAction, CSIDH512Timer, CSIDH512_N, CSIDH512_CURVE_BYTES, CSIDH512_SCALAR_BYTES
from thriso.meter import METER
from thriso.registry import SCHEMES, ORDER, parties_for, signer_set
from thriso.sharing import LiftedShamir, SubgroupShamir

SHARE_NS = [2, 3, 4, 5, 8, 9, 10, 16, 27, 28, 37, 38, 64, 100, 256, 1000, 1369, 1370, 4096,
            10 ** 4, 50653, 50654, 10 ** 5, 10 ** 6, 1407181, 1407182, 10 ** 7, 10 ** 8, 10 ** 9]

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")

FAC512 = {3: 1, 37: 1, 1407181: 1, 51593604295295867744293584889: 1,
          31599414504681995853008278745587832204909: 1}


def run_one(key, ga, t, K, msg=b"benchmark message", n=None, S=None, rg=112, rs=71):
    cls = SCHEMES[key]
    if n is None:
        n, t = parties_for(key, t)
    if S is None:
        S = signer_set(key, n, t)
    METER.reset()
    sch = cls(ga, n=n, t=t, K=K, zk_general=rg, zk_special=rs)
    t0 = time.perf_counter()
    sch.keygen()
    keygen_s = time.perf_counter() - t0
    t0 = time.perf_counter()
    sch.setup(S)
    setup_s = time.perf_counter() - t0
    t0 = time.perf_counter()
    sig = sch.sign(msg, S)
    sign_s = time.perf_counter() - t0
    t0 = time.perf_counter()
    ok = sig is not None and sch.verify(msg, sig)
    verify_s = time.perf_counter() - t0
    summ = METER.summary(ga.curve_bytes, ga.scalar_bytes)
    empty = {"group_actions": 0, "latency_group_actions": 0, "max_party_group_actions": 0, "bytes": 0, "rounds": 0}
    return {
        "scheme": key,
        "label": cls.label,
        "n": n,
        "t": t,
        "signers": len(S),
        "K": K,
        "T": sch.T,
        "valid": ok,
        "keygen": summ.get("keygen", empty),
        "setup": summ.get("setup", empty),
        "sign": summ.get("sign", empty),
        "verify": summ.get("verify", empty),
        "seconds": {"keygen": keygen_s, "setup": setup_s, "sign": sign_s, "verify": verify_s},
        "nonzero_challenges": sum(1 for c in sig[0] if c != 0) if sig is not None else None,
        "signature_bytes": sch.sign_signature_bytes(),
        "pk_bytes": sch.pk_bytes(),
        "secret_bits_per_party": sch.secret_bits_per_party(),
        "active": cls.active,
        "robust": cls.robust,
        "assumptions": cls.assumptions,
        "model": cls.model,
        "zk_general": rg,
        "zk_special": rs,
        "committee_bits_per_party": sch.committee_bits_per_party() if hasattr(sch, "committee_bits_per_party") else 0,
    }


def sweep_t(ts, K, keys):
    ga = MockAction()
    out = []
    for t in ts:
        for key in keys:
            if key == "grassplus" and t > 12:
                continue
            r = run_one(key, ga, t, K)
            out.append(r)
            print("sweep_t", key, t, K, r["valid"], r["sign"]["group_actions"], round(r["seconds"]["sign"], 2), flush=True)
    return out


def sweep_K(t, Ks, keys):
    ga = MockAction()
    out = []
    for K in Ks:
        for key in keys:
            r = run_one(key, ga, t, K)
            out.append(r)
            print("sweep_K", key, t, K, r["valid"], r["sign"]["group_actions"], flush=True)
    return out


def sweep_R(t, K, keys, pairs):
    ga = MockAction()
    out = []
    for rg, rs in pairs:
        for key in keys:
            r = run_one(key, ga, t, K, rg=rg, rs=rs)
            out.append(r)
            print("sweep_R", key, t, K, rg, r["valid"], r["sign"]["group_actions"], flush=True)
    return out


def grass_rss(ts):
    ga = MockAction()
    out = []
    for t in ts:
        n = 2 * t - 1
        cls = SCHEMES["grass"]
        METER.reset()
        sch = cls(ga, n=n, t=t, K=16)
        sch.keygen = _dealer_grass(sch)
        sch.keygen()
        S = list(range(1, t + 1))
        sig = sch.sign(b"m", S)
        ok = sig is not None and sch.verify(b"m", sig)
        s = METER.summary(ga.curve_bytes, ga.scalar_bytes)["sign"]
        out.append({"t": t, "n": n, "shards": sch.M, "valid": ok, "sign": s,
                    "secret_bits_per_party": sch.secret_bits_per_party()})
        print("grass_rss", t, n, sch.M, ok, s["group_actions"], flush=True)
    return out


def _dealer_grass(sch):
    def keygen():
        METER.set_phase("keygen")
        ga = sch.ga
        sch.shard_values = []
        sch.inter = []
        sch.pk = []
        for l in range(sch.K):
            vals = [ga.rand() for _ in sch.shards]
            chain = [ga.E0]
            acc = 0
            for v in vals:
                acc = (acc + v) % ga.N
                chain.append(ga.act(ga.E0, acc))
            sch.shard_values.append(vals)
            sch.inter.append(chain)
            sch.pk.append(chain[-1])
    return keygen


def real_validation(ts, K, keys):
    ga = ToyCSIDHAction()
    out = []
    for t in ts:
        for key in keys:
            r = run_one(key, ga, t, K)
            r["backend"] = "toy-csidh"
            out.append(r)
            print("real", key, t, r["valid"], r["sign"]["group_actions"], round(r["seconds"]["sign"], 1), flush=True)
    return out


def cheat_detection(t, K, trials):
    ga = MockAction()
    rng = random.Random(2026)
    n = 2 * t - 1
    cls = SCHEMES["ours"]
    sch = cls(ga, n=n, t=t, K=K)
    sch.keygen()
    S = list(range(1, t + 1))
    sch.setup(S)
    out = {"t": t, "K": K, "trials": trials}
    honest = 0
    for _ in range(trials):
        sig = sch.sign(b"honest", S)
        honest += int(sig is not None and sch.verify(b"honest", sig))
    out["honest_valid"] = honest
    for typ in ["setup", "chain", "response"]:
        stats = {"abort": 0, "invalid": 0, "valid": 0, "by_position": {}}
        for _ in range(trials):
            pos = rng.randrange(t)
            cheat = {"party": S[pos], "type": typ, "position": rng.randrange(sch.T)}
            if typ == "setup":
                res = "abort" if sch.setup(S, cheat) is None else "valid"
                sch.setup(S)
            else:
                sig = sch.sign(b"cheat", S, cheat)
                if sig is None:
                    res = "abort"
                elif sch.verify(b"cheat", sig):
                    res = "valid"
                else:
                    res = "invalid"
            stats[res] += 1
            key = str(pos + 1)
            stats["by_position"].setdefault(key, {"abort": 0, "invalid": 0, "valid": 0})
            stats["by_position"][key][res] += 1
        out[typ] = stats
        print("cheat", typ, stats["abort"], stats["invalid"], stats["valid"], flush=True)
    return out


def csidh512_timing(samples):
    tm = CSIDH512Timer()
    vals = []
    for _ in range(samples):
        t0 = time.perf_counter()
        tm.one()
        vals.append(time.perf_counter() - t0)
    return {"samples": samples, "mean": statistics.mean(vals), "stdev": statistics.pstdev(vals), "values": vals}


def share_sizes(ns):
    out = []
    for n in ns:
        L = LiftedShamir(CSIDH512_N, FAC512, n, 2)
        row = {"n": n, "lifted_bits": L.share_bits(), "full_bits": math.log2(CSIDH512_N),
               "degrees": [[p, d] for (p, e, d, q, R, u) in L.components]}
        G = SubgroupShamir(CSIDH512_N, FAC512, n, 2)
        row["subgroup_bits"] = G.share_bits()
        row["subgroup_loss_bits"] = G.security_loss_bits()
        out.append(row)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--skip-real", action="store_true")
    args = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    keys = ORDER
    if args.quick:
        ts = [2, 3, 4]
        Ks = [1, 16]
        real_ts = [2]
        samples = 3
        rss_ts = [2, 3]
    else:
        ts = [2, 3, 4, 5, 6, 8, 10, 12, 16]
        Ks = [1, 16, 256]
        real_ts = [2, 3]
        samples = 20
        rss_ts = [2, 3, 4, 5]
    data = {"sizes": {"curve_bytes": CSIDH512_CURVE_BYTES, "scalar_bytes": CSIDH512_SCALAR_BYTES}}

    def save():
        with open(os.path.join(RESULTS, "bench.json"), "w") as fh:
            json.dump(data, fh, indent=1)

    data["csidh512_ga_seconds"] = csidh512_timing(samples)
    data["share_sizes"] = share_sizes(SHARE_NS)
    save()
    data["sweep_t"] = sweep_t(ts, 16, keys)
    data["sweep_K"] = sweep_K(4, Ks, keys)
    data["sweep_R"] = sweep_R(4, 16, keys, [(112, 71), (256, 162)])
    data["grass_rss"] = grass_rss(rss_ts)
    save()
    data["cheat"] = cheat_detection(4, 16, 60 if args.quick else 600)
    save()
    if not args.skip_real:
        data["real"] = real_validation(real_ts, 16, keys)
        save()


if __name__ == "__main__":
    main()
