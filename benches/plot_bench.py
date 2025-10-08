import json
import matplotlib.pyplot as plt
import numpy as np


def r2(y, y_pred):
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot
    return r2


def main():
    filename = "result.json"
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
    print(time_nano_second)

    les_x = range(1, len(time_nano_second) + 1)

    plt.plot(les_x, time_nano_second, "X--", label="raw data")

    coef = np.polyfit(les_x, time_nano_second, 1)
    poly1d_fn = np.poly1d(coef)
    time_nano_second_pred = poly1d_fn(les_x)
    # poly1d_fn is now a function which takes in x and returns an estimate for y
    plt.plot(
        les_x,
        time_nano_second_pred,
        "-",
        label=f"linear regression (R2= {r2(time_nano_second, time_nano_second_pred)})",
    )
    plt.xlabel("number of hash function")
    plt.ylabel("time (ns)")
    plt.ymin = 0
    plt.legend(loc="best")
    plt.title(
        "Time required to index 100000 bases with regars to the number of hash functions"
    )
    plt.show()


if __name__ == "__main__":
    main()
