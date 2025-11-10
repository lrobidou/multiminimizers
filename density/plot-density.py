import json
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from math import log
import matplotlib as mpl

mpl.rcParams["figure.dpi"] = 400
# show more dataframe rows
# pd.set_option("display.min_rows", 50)
# pd.set_option("display.max_rows", 100)


def g(w, k):
    """Lower bound (https://doi.org/10.1101/2024.09.06.611668)"""
    return 1 / (w + k) * ((w + k + (w - 1)) // w)


def read(file):
    with open(file, "r") as f:
        data = json.load(f)
    return pd.json_normalize(data)


def plot(data):
    if isinstance(data, str):
        data = read(data)
    df = data
    s = df["sigma"].unique()[0]
    ws = df.w.unique()

    for w in ws:
        plt.axhline(y=(1) / (w), color="black", linewidth=0.5)
        plt.axhline(y=(2) / (w + 1), color="black", linewidth=0.5)
        ks = range(df.k.min(), df.k.max())
        plt.plot(ks, [g(w, k) for k in ks], color="red", linewidth=0.5)

    sns.lineplot(
        x="k",
        y="density",
        hue="minimizer_name",
        size="w",
        sizes=(1, 2),
        data=df,
        legend="full",
        marker=".",
        dashes=False,
    )
    plt.title(f"Minimizer density on random text (alphabet size σ={s}; length=5M)")
    plt.xlabel("Kmer length k")
    plt.ylabel("Density")
    plt.yscale("log", base=2)
    plt.yticks(
        [2 / (w + 1) for w in ws] + [1 / w for w in ws],
        [f"{2/(w+1):.3f}" for w in ws] + [f"{1/w:.3f}" for w in ws],
    )
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.savefig(f"density_{s}.svg", bbox_inches="tight")
    plt.show()


plot("density_4.json")
