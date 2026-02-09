import json
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from math import log
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import argparse
from pathlib import Path
from matplotlib.ticker import MaxNLocator

plt.switch_backend("agg")


TOOLS_PALETTE = {
    "RandomMinimizer": "tab:blue",
    "ModMinimizer": "tab:green",
    "DoubleDecycling": "tab:red",
    "Miniception": "tab:orange",
    "GreedyMini": "tab:purple",
}

PALETTE = [
    "#4E79A7",
    "#76B7B2",
    "#B07AA1",
    "#FF9DA7",
    "#E15759",
    "#F28E2B",
    "#EDC948",
    "#9C755F",
]

LABEL_ORDER = [
    "RandomMinimizer",
    "ModMinimizer",
    "Miniception",
    "DoubleDecycling",
    "GreedyMini",
    "Lower",
    "Multiminimizers",
    "MOCMM",
]

OCMM_DIR = Path("../multimodmini")


def parse_tuple(s):
    """
    "(1.0,2.0)" -> (1.0, 2.0)"
    """
    s = s.strip()[1:-1]
    x_str, y_str = s.split(",")
    x = float(x_str)
    y = float(y_str)
    return (x, y)


def build_ocmm():
    import subprocess
    import sys

    # Path to the Rust ocmm project
    project_path = OCMM_DIR

    # Compile the Rust project
    result = subprocess.run(
        ["cargo", "build", "--release"],
        cwd=project_path,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print("Compilation failed:")
        print(result.stderr)
        sys.exit(1)


def ocmm(minimizer_size, w, nb_hash):
    import subprocess

    # Path to the compiled binary
    binary_path = project_path / "target" / "release" / project_path.name

    # Execute the binary
    run_result = subprocess.run(
        [
            str(binary_path),
            "--minimizer-size",
            str(minimizer_size),
            "--window",
            str(w),
            "--hashes",
            str(nb_hash),
            "--verify",
            "--size-read",
            str(50000),
        ],
        capture_output=True,
        text=True,
    )

    # print(run_result.stderr)
    assert not run_result.stderr
    res = parse_tuple(run_result.stdout.strip())
    return res[0]


def multimini_load(filename):
    with open(filename, "r") as fichier:
        content_all_k = json.load(fichier)
        ns = content_all_k["ns"]
        content_all_k = content_all_k["data"]
    ws = []
    densities = [[] for _ in range(len(ns))]
    minimizer_sizes = []
    ks = content_all_k.keys()
    ks = [int(k) for k in ks]
    ks.sort()
    for k in ks:
        # get values from the json
        content = content_all_k[str(k)]
        size_read = content["size_read"]
        k = content["k"]
        m = content["m"]
        nb_superkmer = content["nb_superkmer"]

        minimizer_sizes.append(m)

        # compute the density
        for i in range(len(densities)):
            densities[i].append(nb_superkmer[i] / (size_read - m + 1))
        # compute w
        ws.append(k - m + 1)

    return ns, ws, densities, minimizer_sizes


def is_sorted(l):
    return all(a <= b for a, b in zip(l, l[1:]))


def multimini(filename):
    ns, ws, densities, minimizer_sizes = multimini_load(filename)
    for i in range(len(ws)):
        if ws[i] != ws[0]:
            print("internal error")
            exit()
    assert is_sorted(minimizer_sizes)
    return ns, densities, minimizer_sizes


def simple_lower_bound(w, k):
    """Lower bound for forward schemes (theorem 1 from https://doi.org/10.1101/2024.09.06.611668)"""
    return ((w + k + w - 1) // w) / (w + k)


def lower_bound(w, k):
    """Improved lower bound for forward schemes (theorem 2 from https://doi.org/10.1101/2024.09.06.611668)"""
    k_ = ((k - 1 + w - 1) // w) * w + 1
    return max(simple_lower_bound(w, k), simple_lower_bound(w, k_))


def read(file):
    with open(file, "r") as f:
        data = json.load(f)
    return pd.json_normalize(data)


def plot(data, plot_mocmm: bool, plot_only_small_values: bool):
    if isinstance(data, str):
        data = read(data)
    df = data
    s = df["sigma"].unique()[0]
    ws = df.w.unique()

    from matplotlib import rc

    # rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
    # rc("text", usetex=True)
    # fontsize = 30

    for w in ws:
        df_filtered = df[df["w"] == w].drop(columns=["w"])
        plt.axhline(y=(1) / (w), color="black", linewidth=0.5)
        plt.axhline(y=(2) / (w + 1), color="black", linewidth=0.5)
        ks = list(range(df_filtered.k.min(), df_filtered.k.max() + 1))

        sns.lineplot(
            x="k",
            y="density",
            hue="minimizer_name",
            sizes=(1, 2),
            data=df_filtered,
            legend="full",
            marker=".",
            dashes=False,
            palette=TOOLS_PALETTE,
        )

        plt.plot(
            ks,
            [lower_bound(w, k) for k in ks],
            linewidth=2,
            color="black",
            label="Lower bound for forward schemes",
        )

        ns, densities, minimizer_sizes = multimini(f"../data_fixed_w_{w}.json")
        for nb_hash, nb_hash_data, color in zip(ns, densities, PALETTE):
            if not plot_mocmm:
                if nb_hash > 1:
                    plt.plot(
                        minimizer_sizes,
                        nb_hash_data,
                        "--",
                        color=color,
                        label=f"Multiminimizers ({nb_hash})",
                    )
            else:
                ocmms = []

                run_ocmm_partial = partial(ocmm, w=w, nb_hash=nb_hash)

                with ProcessPoolExecutor(79) as executor:
                    ocmms = list(executor.map(run_ocmm_partial, ks))
                plt.plot(
                    ks,
                    ocmms,
                    "--",
                    color=color,
                    label=f"MOCMM ({nb_hash})",
                )

        if plot_only_small_values:
            plt.xlim(ks[0], 31)

        plt.title(f"Minimizer density on random text ($w={w}$)")
        plt.xlabel("Minimizer length $m$")
        plt.ylabel("Density")
        plt.yticks(
            [2 / (w + 1)] + [1 / w],
            [f"{2 / (w + 1):.3f}"] + [f"{1 / w:.3f}"],
        )

        handles, labels = plt.gca().get_legend_handles_labels()
        pairs = list(zip(labels, handles))
        f = lambda p: (LABEL_ORDER.index(p[0].split(" ")[0]), len(p[0]), p[0])
        pairs.sort(key=f)
        labels, handles = zip(*pairs)
        plt.legend(
            handles,
            labels,
            bbox_to_anchor=(0.5, -0.15),
            loc="upper center",
            ncol=2,
        )

        ax = plt.gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.savefig(f"density_{s}_w_{w}.pdf", bbox_inches="tight", dpi=300)
        plt.clf()
        plt.close()
        sns.reset_defaults()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot_mocmm", default=False, action="store_true")
    parser.add_argument("--plot_only_small_values", default=False, action="store_true")
    args = parser.parse_args()

    plot_mocmm = args.plot_mocmm
    plot_only_small_values = args.plot_only_small_values

    if plot_mocmm:
        build_ocmm()

    plot("density_4.json", plot_mocmm, plot_only_small_values)


if __name__ == "__main__":
    main()
