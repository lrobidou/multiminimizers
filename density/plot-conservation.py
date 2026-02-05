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

PALETTE = [
    "#4E79A7",
    "#76B7B2",
    "#B07AA1",
    "#E15759",
    "#F28E2B",
    "#EDC948",
]


def load_conservation(filename):
    with open(filename, "r") as f:
        content = json.load(f)
        m = content["m"]
        w = content["w"]
        error_rates = content["error_rates"]
        data = content["data"]
    ns = sorted(list(map(int, data.keys())))
    conservations = [data[str(n)] for n in ns]
    return m, w, error_rates, ns, conservations


def plot(filename):
    m, w, error_rates, ns, conservations = load_conservation(filename)

    # from matplotlib import rc
    # rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
    # rc("text", usetex=True)

    for nb_hash, conservation, color in zip(ns, conservations, PALETTE):
        plt.plot(
            error_rates,
            conservation,
            # "--",
            color=color,
            label=f"N={nb_hash}",
        )

    plt.xlim(left=0)
    plt.ylim((0, 1))
    plt.grid()
    plt.title(f"Conservation of multiminimizers on random text ($m={m}, w={w}$)")
    plt.xlabel("Error rate")
    plt.ylabel("Conservation (Jaccard similarity)")
    plt.legend()

    # plt.tight_layout()
    plt.savefig(f"conservation.pdf", bbox_inches="tight", dpi=300)
    plt.clf()
    plt.close()


if __name__ == "__main__":
    plot("../data_conservation.json")
