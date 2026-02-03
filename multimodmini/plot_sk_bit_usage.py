import json
import matplotlib.pyplot as plt
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import argparse
from pathlib import Path
from matplotlib import rc
from matplotlib.ticker import MaxNLocator

plt.switch_backend("agg")


PALETTE = [
    "#4E79A7",
    "#76B7B2",
    "#B07AA1",
    "#FF9DA7",
    "#E15759",
    "#F28E2B",
    "#EDC948",
    "#59A14F",
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


def ocmm(minimizer_size, w, nb_hash):
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
    stdout = run_result.stdout.strip()
    # print(stdout)
    res = parse_tuple(stdout)
    return res


def lower_bound(w, k):
    """Lower bound for forward schemes (https://doi.org/10.1101/2024.09.06.611668)"""
    return 1 / (w + k) * ((w + k + (w - 1)) // w)


def read(file):
    with open(file, "r") as f:
        data = json.load(f)
    return pd.json_normalize(data)


def plot(ks):
    m = 21

    rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
    rc("text", usetex=True)
    fontsize = 30

    plt.rcParams["figure.figsize"] = (12, 6)

    for i, k in enumerate(ks):
        nbs_hashes = list(range(1, 17))

        assert k >= m
        w = k - m + 1

        run_ocmm_partial = partial(ocmm, m, w)

        with ProcessPoolExecutor(79) as executor:
            ocmms = list(executor.map(run_ocmm_partial, nbs_hashes))
        ocmm_super_kmer_bit_usage = [x[1] for x in ocmms]

        plt.plot(
            nbs_hashes,
            ocmm_super_kmer_bit_usage,
            "--",
            color=PALETTE[i],
            label=f"k = {k}",
        )

    # plt.title(f"Minimizer density on random text ($w={w}$)", fontsize=fontsize)
    plt.xlabel("Number of hash functions", fontsize=fontsize)
    plt.ylabel("Bits / $k$-mers", fontsize=fontsize)

    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend(loc="upper right", ncols=2, fontsize=fontsize)
    plt.tick_params(axis="both", which="major", labelsize=fontsize)

    plt.tight_layout()
    plt.savefig("bit_usage_superkmer_mocmm.pdf")
    plt.clf()
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument("--plot_only_small_values", default=False, action="store_true")
    parser.add_argument("-k", nargs="+", required=True)
    args = parser.parse_args()

    # plot_mocmm = args.plot_mocmm
    # plot_only_small_values = args.plot_only_small_values
    ks = args.k
    ks = [int(k) for k in ks]

    plot(ks)


if __name__ == "__main__":
    main()
