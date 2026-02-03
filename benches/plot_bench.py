import json
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib import rc
from matplotlib.ticker import MaxNLocator


def r2(y, y_pred):
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot
    return r2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    args = parser.parse_args()

    filename = args.filename

    with open(filename, "r") as fichier:
        benches = []
        index = 1
        for ligne in fichier:
            bench = json.loads(ligne)
            bench_id = bench.get("id")

            if (bench_id is not None) and bench_id.startswith("Index/"):
                expected_id = f"Index/{index} hash function"
                assert bench_id == expected_id
                benches.append(bench)
                index += 1

    time_nano_second = []
    for bench in benches:
        assert bench["mean"]["unit"] == "ns"
        time_nano_second.append(bench["mean"]["estimate"])

    les_x = range(1, len(time_nano_second) + 1)

    rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
    rc("text", usetex=True)
    fontsize = 30

    plt.rcParams["figure.figsize"] = (12, 6)
    plt.plot(les_x, time_nano_second, "X--", label="raw data")

    coef = np.polyfit(les_x, time_nano_second, 1)
    poly1d_fn = np.poly1d(coef)
    time_nano_second_pred = poly1d_fn(les_x)
    # poly1d_fn is now a function which takes in x and returns an estimate for y
    plt.plot(
        les_x,
        time_nano_second_pred,
        "-",
        label=f"linear regression (R2= {round(r2(time_nano_second, time_nano_second_pred), 5)})",
    )
    plt.xlabel("number of hash function", fontsize=fontsize)
    plt.ylabel("time (ns)", fontsize=fontsize)
    plt.ymin = 0
    # plt.legend(loc="best", fontsize=fontsize)
    plt.title(
        "Time required to iterate over the super-k-mers\nwith regard to the number of hash functions",
        fontsize=fontsize,
    )
    # plt.show()
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend(loc="best", fontsize=fontsize)
    plt.tick_params(axis="both", which="major", labelsize=fontsize)

    plt.tight_layout()
    plt.savefig("time_vs_N.pdf")


if __name__ == "__main__":
    main()
