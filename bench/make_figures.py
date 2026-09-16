import argparse
import json
import math
import os
import shutil
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from thriso.groupaction import CSIDH512_N
from thriso.sharing import LiftedShamir, SubgroupShamir

RESULTS = os.path.join(ROOT, "results")
FIGDIR = os.path.join(RESULTS, "figures")

FAC512 = {3: 1, 37: 1, 1407181: 1, 51593604295295867744293584889: 1,
          31599414504681995853008278745587832204909: 1}

ORDER = ["dm20", "dems24", "sashimi", "csishark", "cm22", "thresher", "grass", "grassplus", "ours"]
LABEL = {"dm20": "DM20", "dems24": "DEMS24", "sashimi": "CS20", "csishark": "ABCP23a", "cm22": "CM22",
         "thresher": "ABCP23b", "grass": "BBMP24", "grassplus": "BBDMP25", "ours": "Albacore"}
COLOR = {"dm20": "#7f7f7f", "dems24": "#b5b04a", "sashimi": "#1f77b4", "csishark": "#17becf", "cm22": "#9467bd",
         "thresher": "#8c564b", "grass": "#2ca02c", "grassplus": "#ff7f0e", "ours": "#d62728"}
MARKER = {"dm20": "o", "dems24": "v", "sashimi": "s", "csishark": "x", "cm22": "D", "thresher": "^",
          "grass": "P", "grassplus": "h", "ours": "*"}
LINE = {"dm20": ":", "dems24": ":", "sashimi": "-", "csishark": "--", "cm22": "-", "thresher": "-",
        "grass": "-", "grassplus": "-", "ours": "-"}


def setup_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "ps.fonttype": 42,
        "pdf.fonttype": 42,
        "axes.grid": True,
        "grid.linestyle": ":",
        "grid.color": "0.75",
        "grid.linewidth": 0.6,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.2,
        "lines.markersize": 4.5,
        "savefig.dpi": 300,
    })


def save(fig, name, paper_dir, tight=True):
    os.makedirs(FIGDIR, exist_ok=True)
    eps = os.path.join(FIGDIR, name + ".eps")
    pdf = os.path.join(FIGDIR, name + ".pdf")
    if tight:
        fig.savefig(eps, format="eps", bbox_inches="tight", pad_inches=0.03)
    else:
        fig.savefig(eps, format="eps")
    plt.close(fig)
    subprocess.run(["epstopdf", eps, "--outfile=" + pdf], check=True)
    if paper_dir:
        os.makedirs(paper_dir, exist_ok=True)
        shutil.copy(pdf, os.path.join(paper_dir, name + ".pdf"))
    print("figure", name)


def index(rows):
    return {(r["scheme"], r["t"]): r for r in rows}


def style_kw(key):
    kw = dict(color=COLOR[key], marker=MARKER[key], linestyle=LINE[key],
              linewidth=2.0 if key == "ours" else 1.2,
              markersize=8 if key == "ours" else 4.5, label=LABEL[key],
              zorder=4 if key == "ours" else 3)
    if key in ("dm20", "dems24", "grass"):
        kw.update(markerfacecolor="white", markeredgewidth=1.0, markersize=5, zorder=6)
    return kw


def series(rows, key, field):
    pts = sorted((r["t"], field(r)) for r in rows if r["scheme"] == key)
    return [p[0] for p in pts], [p[1] for p in pts]


def fig_sign(data, paper_dir):
    rows = data["sweep_t"]
    tau = data["csidh512_ga_seconds"]["mean"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.75))
    for key in ORDER:
        x, y = series(rows, key, lambda r: r["sign"]["group_actions"])
        axes[0].plot(x, y, **style_kw(key))
        x, y = series(rows, key, lambda r: r["sign"]["latency_group_actions"] * tau)
        axes[1].plot(x, y, **style_kw(key))
    rss = data.get("grass_rss", [])
    if rss:
        axes[0].plot([r["t"] for r in rss], [r["sign"]["group_actions"] for r in rss], color=COLOR["grass"],
                     marker=MARKER["grass"], linestyle="--", linewidth=1.0, label="BBMP24, $n=2t-1$")
        axes[1].plot([r["t"] for r in rss], [r["sign"]["latency_group_actions"] * tau for r in rss],
                     color=COLOR["grass"], marker=MARKER["grass"], linestyle="--", linewidth=1.0)
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")
    axes[0].set_xlabel("number of signers $t$")
    axes[1].set_xlabel("number of signers $t$")
    axes[0].set_ylabel("group actions per signature")
    axes[1].set_ylabel("latency on CSIDH-512 (s)")
    axes[0].set_title("(a) total group actions")
    axes[1].set_title("(b) latency with parallel parties")
    for ax in axes:
        ax.set_xticks([2, 4, 6, 8, 10, 12, 14, 16])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.12), frameon=False,
               columnspacing=1.2, handlelength=2.2)
    fig.tight_layout()
    save(fig, "fig_sign", paper_dir)


def fig_comm_state(data, paper_dir):
    rows = data["sweep_t"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.75))
    for key in ORDER:
        x, y = series(rows, key, lambda r: r["sign"]["bytes"] / 1000.0)
        axes[0].plot(x, y, **style_kw(key))
        x, y = series(rows, key, lambda r: (r["secret_bits_per_party"] + r.get("committee_bits_per_party", 0)) / 8000.0)
        axes[1].plot(x, y, **style_kw(key))
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")
    axes[0].set_xlabel("number of signers $t$")
    axes[1].set_xlabel("number of signers $t$")
    axes[0].set_ylabel("signing communication (kB)")
    axes[1].set_ylabel("secret state per party (kB)")
    axes[0].set_title("(a) communication per signature")
    axes[1].set_title("(b) persistent and committee state")
    for ax in axes:
        ax.set_xticks([2, 4, 6, 8, 10, 12, 14, 16])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=9, bbox_to_anchor=(0.5, 1.07), frameon=False,
               columnspacing=0.9, handlelength=2.0)
    fig.tight_layout()
    save(fig, "fig_comm_state", paper_dir)


def transition_points(limit):
    ns = {2}
    for p in [3, 37, 1407181]:
        q = p
        while q <= limit:
            ns.add(q)
            ns.add(q + 1)
            q *= p
    ns.add(limit)
    return sorted(n for n in ns if 2 <= n <= limit)


def fig_share(data, paper_dir):
    limit = 10 ** 9
    ns = transition_points(limit)
    lifted = []
    subgroup = []
    for n in ns:
        lifted.append(LiftedShamir(CSIDH512_N, FAC512, n, 2).share_bits())
        subgroup.append(SubgroupShamir(CSIDH512_N, FAC512, n, 2).share_bits())
    full = math.log2(CSIDH512_N)
    fig, ax = plt.subplots(1, 1, figsize=(4.6, 2.6))
    ax.step(ns, lifted, where="post", color=COLOR["ours"], linewidth=1.8, label="lifted share size (Albacore)")
    ax.step(ns, subgroup, where="post", color=COLOR["dm20"], linewidth=1.4, linestyle="--",
            label="subgroup key space (DM20)")
    ax.axhline(full, color="0.3", linewidth=0.9, linestyle=":", label=r"$\log N$")
    meas = data.get("share_sizes", [])
    if meas:
        ax.plot([r["n"] for r in meas], [r["lifted_bits"] for r in meas], linestyle="none", marker="o",
                markersize=3.2, markerfacecolor="white", markeredgecolor=COLOR["ours"], label="benchmark points")
    for p, name in [(3, "3"), (37, "37"), (1407181, "1407181")]:
        ax.axvline(p, color="0.6", linewidth=0.6, linestyle="-.")
        ax.text(p * 1.15, 334, "$n=" + name + "$", fontsize=7, va="center", ha="left", color="0.35")
    ax.set_xscale("log")
    ax.set_xlim(1.5, limit)
    ax.set_ylim(222, 340)
    ax.set_xlabel("number of parties $n$")
    ax.set_ylabel("bits")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, columnspacing=1.2)
    fig.tight_layout()
    save(fig, "fig_share", paper_dir)


def grouped_bars(ax, groups, keys, values, colors, labels, hatches=None):
    width = 0.8 / len(groups)
    x = np.arange(len(keys))
    for gi, g in enumerate(groups):
        vals = [values[(k, g)] for k in keys]
        ax.bar(x + (gi - (len(groups) - 1) / 2) * width, vals, width, color=colors[gi], edgecolor="black",
               linewidth=0.4, label=labels[gi], hatch=None if hatches is None else hatches[gi], zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([LABEL[k] for k in keys], rotation=35, ha="right")


def fig_params(data, paper_dir):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
    Ks = sorted({r["K"] for r in data["sweep_K"]})
    vals = {(r["scheme"], r["K"]): r["sign"]["group_actions"] for r in data["sweep_K"]}
    grouped_bars(axes[0], Ks, ORDER, vals, ["#c6dbef", "#6baed6", "#08519c"][:len(Ks)],
                 ["$K=%d$, $T=%d$" % (K, next(r["T"] for r in data["sweep_K"] if r["K"] == K)) for K in Ks])
    axes[0].set_yscale("log")
    axes[0].set_ylabel("group actions per signature")
    axes[0].set_title("(a) public curves $K$ at $t=4$")
    axes[0].legend(loc="upper left", ncol=3, frameon=True, framealpha=1.0, edgecolor="0.8", fontsize=6.8)
    axes[0].set_ylim(10, 5e7)
    Rs = sorted({r["zk_general"] for r in data["sweep_R"]})
    valsR = {(r["scheme"], r["zk_general"]): r["sign"]["group_actions"] for r in data["sweep_R"]}
    labelsR = []
    for R in Rs:
        rs = next(r["zk_special"] for r in data["sweep_R"] if r["zk_general"] == R)
        labelsR.append("$R_G=%d$, $R_S=%d$" % (R, rs))
    grouped_bars(axes[1], Rs, ORDER, valsR, ["#fdd0a2", "#d94801"], labelsR)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("group actions per signature")
    axes[1].set_title("(b) proof repetitions at $t=4$, $K=16$")
    axes[1].legend(loc="upper left", ncol=2, frameon=True, framealpha=1.0, edgecolor="0.8", fontsize=6.8)
    axes[1].set_ylim(10, 5e7)
    fig.tight_layout()
    save(fig, "fig_params", paper_dir)


def fig_amort(data, paper_dir):
    idx = index(data["sweep_t"])
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    Q = np.unique(np.round(np.logspace(0, 5, 200)).astype(int))
    for ax, t in zip(axes, [4, 8]):
        for key in ORDER:
            r = idx.get((key, t))
            if r is None:
                continue
            setup = r["setup"]["group_actions"]
            sign = r["sign"]["group_actions"]
            kw = style_kw(key)
            kw["markevery"] = 25
            ax.plot(Q, (setup + Q * sign) / Q, **kw)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("signatures $Q$ per committee")
        ax.set_ylabel("group actions per signature")
        ax.set_title("(%s) $t=%d$, $K=16$, setup included" % ("a" if t == 4 else "b", t))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=9, bbox_to_anchor=(0.5, 1.07), frameon=False,
               columnspacing=0.9, handlelength=2.0)
    fig.tight_layout()
    save(fig, "fig_amort", paper_dir)


def fig_real(data, paper_dir):
    rows = data.get("real")
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1.15, 1]})
    xs = np.array([r["sign"]["group_actions"] for r in rows], dtype=float)
    ys = np.array([r["seconds"]["sign"] for r in rows], dtype=float)
    slope = float((xs * ys).sum() / (xs * xs).sum())
    grid = np.logspace(math.log10(xs.min() * 0.7), math.log10(xs.max() * 1.4), 50)
    axes[0].plot(grid, slope * grid, color="0.4", linewidth=0.9, linestyle="--",
                 label="fit: %.2f ms per group action" % (1000 * slope))
    for key in ORDER:
        pts = [r for r in rows if r["scheme"] == key]
        if not pts:
            continue
        axes[0].plot([r["sign"]["group_actions"] for r in pts], [r["seconds"]["sign"] for r in pts],
                     linestyle="none", marker=MARKER[key], color=COLOR[key],
                     markersize=8 if key == "ours" else 5, label=LABEL[key])
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("counted group actions")
    axes[0].set_ylabel("measured signing time (s)")
    axes[0].set_title("(a) real isogeny backend, $t\\in\\{2,3\\}$")
    axes[0].legend(loc="upper left", ncol=2, frameon=True, framealpha=1.0, edgecolor="0.8", fontsize=6.3)
    ts = sorted({r["t"] for r in rows})
    width = 0.8 / len(ts)
    x = np.arange(len(ORDER))
    for gi, t in enumerate(ts):
        vals = []
        for key in ORDER:
            r = next((r for r in rows if r["scheme"] == key and r["t"] == t), None)
            vals.append(1000 * r["seconds"]["sign"] / r["sign"]["group_actions"] if r else 0)
        axes[1].bar(x + (gi - (len(ts) - 1) / 2) * width, vals, width, color=["#9ecae1", "#3182bd"][gi % 2],
                    edgecolor="black", linewidth=0.4, label="$t=%d$" % t, zorder=3)
    axes[1].axhline(1000 * slope, color="0.4", linewidth=0.9, linestyle="--")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([LABEL[k] for k in ORDER], rotation=35, ha="right")
    axes[1].set_ylabel("ms per group action")
    axes[1].set_title("(b) time per counted group action")
    axes[1].legend(loc="upper right", frameon=True, framealpha=1.0, edgecolor="0.8")
    fig.tight_layout()
    save(fig, "fig_real", paper_dir)


def fig_3d(data, paper_dir):
    idx = index(data["sweep_t"])
    ts = sorted({r["t"] for r in data["sweep_t"]})
    ref = max(ts)
    keys = sorted(ORDER, key=lambda k: idx[(k, 4)]["sign"]["group_actions"])
    fig = plt.figure(figsize=(5.9, 4.5))
    ax = fig.add_subplot(111, projection="3d")
    dx = 0.5
    dy = 0.5
    for yi in reversed(range(len(keys))):
        key = keys[yi]
        for xi in reversed(range(len(ts))):
            r = idx.get((key, ts[xi]))
            if r is None:
                continue
            z = math.log10(r["sign"]["group_actions"])
            ax.bar3d(xi - dx / 2, yi - dy / 2, 0, dx, dy, z, color=COLOR[key], edgecolor="black", linewidth=0.2,
                     shade=True)
    ax.set_xticks(range(len(ts)))
    ax.set_xticklabels([str(t) for t in ts], fontsize=7)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([LABEL[k] for k in keys], fontsize=7, ha="left", va="center")
    ax.set_zlim(0, 7)
    ax.set_zticks([0, 1, 2, 3, 4, 5, 6, 7])
    ax.set_zticklabels(["$10^{%d}$" % v for v in range(0, 8)], fontsize=7)
    ax.set_xlabel("number of signers $t$", labelpad=2)
    ax.set_zlabel("group actions per signature", labelpad=7, fontsize=8.5)
    ax.view_init(elev=20, azim=-128)
    ax.set_box_aspect((1.3, 1.3, 0.75))
    fig.subplots_adjust(left=0.05, right=1.03, bottom=-0.02, top=1.06)
    save(fig, "fig_3d", paper_dir, tight=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper-dir", default=None)
    args = ap.parse_args()
    setup_style()
    with open(os.path.join(RESULTS, "bench.json")) as fh:
        data = json.load(fh)
    if "sweep_t" in data:
        fig_sign(data, args.paper_dir)
        fig_comm_state(data, args.paper_dir)
        fig_params(data, args.paper_dir)
        fig_amort(data, args.paper_dir)
        fig_3d(data, args.paper_dir)
    fig_share(data, args.paper_dir)
    fig_real(data, args.paper_dir)


if __name__ == "__main__":
    main()
