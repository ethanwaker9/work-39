import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
TABDIR = os.path.join(RESULTS, "tables")

ORDER = ["dm20", "dems24", "sashimi", "csishark", "cm22", "thresher", "grass", "grassplus", "ours"]
LABEL = {"dm20": "DM20", "dems24": "DEMS24", "sashimi": "CS20", "csishark": "ABCP23a", "cm22": "CM22",
         "thresher": "ABCP23b", "grass": "BBMP24", "grassplus": "BBDMP25", "ours": "\\Albacore"}
CITE = {"dm20": "DM20", "dems24": "DEMS24", "sashimi": "CS20", "csishark": "ABCP23a", "cm22": "CM22",
        "thresher": "ABCP23b", "grass": "BBMP24", "grassplus": "BBDMP25"}


def sign_formula(key, t, n, T, K, RS, RG, Tne, M=None):
    if key == "dm20":
        return t * T
    if key == "dems24":
        return t * T + K
    if key in ("sashimi", "csishark"):
        return T * (2 * t + (t * t + 2 * t) * RS + 2 * t * (t - 1) * RG)
    if key == "cm22":
        return T * (3 * t + 3 * t * t * RG)
    if key == "thresher":
        return T * (n + n * n * RG)
    if key == "grass":
        if M is None:
            M = math.comb(n, t - 1)
        return T * M * (t + 1)
    if key == "grassplus":
        return T * (t + t * (RS + (t - 1) * RG) + (t - 1) * RS + (t - 1) * (t - 2) // 2 * RG + t * (t - 1) // 2)
    if key == "ours":
        return t * T + T + (t - 1) * Tne
    raise KeyError(key)


def setup_formula(key, t, K, RG):
    if key == "ours":
        return K * (t + (t - 1) * (t + 2) * RG // 2)
    return 0


def check_formulas(data):
    bad = 0
    total = 0
    groups = [("sweep_t", data.get("sweep_t", [])), ("sweep_K", data.get("sweep_K", [])),
              ("sweep_R", data.get("sweep_R", [])), ("real", data.get("real", []))]
    for name, rows in groups:
        for r in rows:
            total += 1
            f = sign_formula(r["scheme"], r["t"], r["n"], r["T"], r["K"], r["zk_special"], r["zk_general"],
                             r.get("nonzero_challenges") or 0)
            s = setup_formula(r["scheme"], r["t"], r["K"], r["zk_general"])
            ok = f == r["sign"]["group_actions"] and s == r["setup"]["group_actions"] and r["valid"]
            if not ok:
                bad += 1
                print("mismatch", name, r["scheme"], r["t"], r["K"], r["sign"]["group_actions"], f,
                      r["setup"]["group_actions"], s, r["valid"])
    for r in data.get("grass_rss", []):
        total += 1
        f = sign_formula("grass", r["t"], r["n"], 23, 16, 71, 112, 0)
        if f != r["sign"]["group_actions"] or not r["valid"]:
            bad += 1
            print("mismatch grass_rss", r["t"], r["sign"]["group_actions"], f)
    print("formula checks", total - bad, "of", total)
    return total, bad


def fmt_int(v):
    return "{:,}".format(int(v)).replace(",", "\\,")


def fmt_float(v, digits=2):
    if v >= 1000:
        return fmt_int(round(v))
    return ("{:.%df}" % digits).format(v)


def table_main(data, t=4):
    tau = data["csidh512_ga_seconds"]["mean"]
    rows = {r["scheme"]: r for r in data["sweep_t"] if r["t"] == t}
    lines = []
    for key in ORDER:
        r = rows[key]
        f = sign_formula(key, r["t"], r["n"], r["T"], r["K"], r["zk_special"], r["zk_general"],
                         r.get("nonzero_challenges") or 0)
        name = LABEL[key] if key == "ours" else "%s \\cite{%s}" % (LABEL[key], CITE[key])
        state = (r["secret_bits_per_party"] + r.get("committee_bits_per_party", 0)) / 8000.0
        lines.append(" & ".join([
            name,
            "$%d$" % r["n"],
            "$%s$" % fmt_int(r["sign"]["group_actions"]),
            "$%s$" % fmt_int(f),
            "$%s$" % fmt_int(r["sign"]["latency_group_actions"]),
            "$%s$" % fmt_float(r["sign"]["latency_group_actions"] * tau / 60.0, 2),
            "$%s$" % fmt_float(r["sign"]["bytes"] / 1000.0, 1),
            "$%d$" % r["sign"]["rounds"],
            "$%s$" % fmt_float(state, 2),
            "$%s$" % fmt_int(r["setup"]["group_actions"]),
        ]) + "\\\\")
        if key == "grass":
            rss = next((x for x in data.get("grass_rss", []) if x["t"] == t), None)
            if rss is not None:
                f = sign_formula("grass", rss["t"], rss["n"], r["T"], r["K"], r["zk_special"], r["zk_general"], 0)
                lines.append(" & ".join([
                    "%s \\cite{%s}, replicated" % (LABEL[key], CITE[key]),
                    "$%d$" % rss["n"],
                    "$%s$" % fmt_int(rss["sign"]["group_actions"]),
                    "$%s$" % fmt_int(f),
                    "$%s$" % fmt_int(rss["sign"]["latency_group_actions"]),
                    "$%s$" % fmt_float(rss["sign"]["latency_group_actions"] * tau / 60.0, 2),
                    "$%s$" % fmt_float(rss["sign"]["bytes"] / 1000.0, 1),
                    "$%d$" % rss["sign"]["rounds"],
                    "$%s$" % fmt_float(rss["secret_bits_per_party"] / 8000.0, 2),
                    "$0$",
                ]) + "\\\\")
    return "\n".join(lines)


def table_real(data):
    rows = data.get("real", [])
    if not rows:
        return ""
    ts = sorted({r["t"] for r in rows})
    lines = []
    for key in ORDER:
        cells = [LABEL[key] if key == "ours" else "%s \\cite{%s}" % (LABEL[key], CITE[key])]
        for t in ts:
            r = next((x for x in rows if x["scheme"] == key and x["t"] == t), None)
            if r is None:
                cells += ["--", "--", "--"]
                continue
            cells += ["$%s$" % fmt_int(r["sign"]["group_actions"]),
                      "$%s$" % fmt_float(r["seconds"]["sign"], 1),
                      "$%s$" % fmt_float(1000 * r["seconds"]["sign"] / r["sign"]["group_actions"], 2)]
        lines.append(" & ".join(cells) + "\\\\")
    return "\n".join(lines)


def table_cheat(data):
    c = data.get("cheat")
    if not c:
        return ""
    t = c["t"]
    lines = []
    names = {"setup": "committee setup", "chain": "commitment chain", "response": "response chain"}
    for typ in ["setup", "chain", "response"]:
        st = c[typ]
        cells = [names[typ]]
        for pos in range(1, t + 1):
            b = st["by_position"].get(str(pos), {"abort": 0, "invalid": 0, "valid": 0})
            cells.append("$%d/%d/%d$" % (b["abort"], b["invalid"], b["valid"]))
        cells.append("$%d/%d/%d$" % (st["abort"], st["invalid"], st["valid"]))
        lines.append(" & ".join(cells) + "\\\\")
    return "\n".join(lines)


def table_amort(data, ts=(2, 4, 8, 16)):
    idx = {(r["scheme"], r["t"]): r for r in data["sweep_t"]}
    lines = []
    for t in ts:
        o = idx.get(("ours", t))
        if o is None:
            continue
        cells = ["$%d$" % t, "$%s$" % fmt_int(o["setup"]["group_actions"]), "$%s$" % fmt_int(o["sign"]["group_actions"])]
        for key in ["grass", "grassplus", "cm22", "sashimi"]:
            r = idx.get((key, t))
            if r is None:
                cells.append("--")
                continue
            gap = r["sign"]["group_actions"] - o["sign"]["group_actions"]
            if gap <= 0:
                cells.append("never")
            else:
                q = o["setup"]["group_actions"] // gap + 1
                cells.append("$%s$" % fmt_int(q))
        lines.append(" & ".join(cells) + "\\\\")
    return "\n".join(lines)


def main():
    with open(os.path.join(RESULTS, "bench.json")) as fh:
        data = json.load(fh)
    os.makedirs(TABDIR, exist_ok=True)
    total, bad = check_formulas(data)
    out = {"tab_main.tex": table_main(data), "tab_real.tex": table_real(data), "tab_cheat.tex": table_cheat(data),
           "tab_amort.tex": table_amort(data)}
    for name, body in out.items():
        with open(os.path.join(TABDIR, name), "w") as fh:
            fh.write(body + "\n")
        print("table", name)
    with open(os.path.join(TABDIR, "summary.json"), "w") as fh:
        tau = data["csidh512_ga_seconds"]
        json.dump({"formula_checks": total, "formula_mismatches": bad, "csidh512_ga_mean_s": tau["mean"],
                   "csidh512_ga_stdev_s": tau["stdev"]}, fh, indent=1)
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
