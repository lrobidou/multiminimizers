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


def multimini_load(filename):
    with open(filename, "r") as fichier:
        content_all_k = json.load(fichier)
        content_all_k = content_all_k["data"]
    ws = []
    densities = [[] for _ in range(8)]
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

    return ws, densities, minimizer_sizes


def is_sorted(l):
    return all(a <= b for a, b in zip(l, l[1:]))


def multimini(filename):
    ws, densities, minimizer_sizes = multimini_load(filename)
    for i in range(len(ws)):
        if ws[i] != ws[0]:
            print("internal error")
            exit()

    assert is_sorted(minimizer_sizes)
    assert is_sorted(minimizer_sizes)

    return densities, minimizer_sizes


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
        df_filtered = df[df["w"] == w].drop(columns=["w"])
        plt.axhline(y=(1) / (w), color="black", linewidth=0.5)
        plt.axhline(y=(2) / (w + 1), color="black", linewidth=0.5)
        ks = range(df_filtered.k.min(), df_filtered.k.max())
        plt.plot(ks, [g(w, k) for k in ks], color="red", linewidth=0.5)

        sns.lineplot(
            x="k",
            y="density",
            hue="minimizer_name",
            # size="w",
            sizes=(1, 2),
            data=df_filtered,
            legend="full",
            marker=".",
            dashes=False,
        )

        densities, minimizer_sizes = multimini(f"../data_fixed_w_{w}.json")
        for nb_hash, nb_hash_data in enumerate(densities):
            plt.plot(minimizer_sizes, nb_hash_data, label=f"multiminimi {nb_hash + 1}")

        plt.title(
            f"Minimizer density on random text (alphabet size σ={s}; length=5M, w={w})"
        )
        plt.xlabel("Minimizer length m")
        plt.ylabel("Density")
        plt.yscale("log", base=2)
        plt.yticks(
            [2 / (w + 1)] + [1 / w],
            [f"{2 / (w + 1):.3f}"] + [f"{1 / w:.3f}"],
        )
        plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

        plt.savefig(f"density_{s}_w_{w}.pdf", bbox_inches="tight")
        plt.clf()
        plt.close()
        sns.reset_defaults()


plot("density_4.json")
