import matplotlib.pyplot as plt
import json

filenames = ["data_canonical.txt", "data_non_canonical.txt"]


def range1(collection):
    return range(1, len(collection) + 1)


for filename in filenames:
    with open(filename, "r") as fichier:
        # plt.figure()
        content = json.load(fichier)

        size_read = content["size_read"]
        k = content["k"]
        m = content["m"]
        limit_overhead = content["limit_overhead"]
        overhead = content["overhead"]
        limit_avg_superkmer_size = content["limit_avg_superkmer_size"]
        avg_superkmer_size = content["avg_superkmer_size"]
        limit_nb_superkmer = content["limit_nb_superkmer"]
        nb_superkmer = content["nb_superkmer"]

        # compute the density
        density = [x * k / size_read for x in nb_superkmer]

        # Compute percentage differences from the first element
        overhead_diff_pct = [
            (x - overhead[0]) / overhead[0] * 100 + 100 for x in overhead
        ]
        avg_size_diff_pct = [
            (x - avg_superkmer_size[0]) / avg_superkmer_size[0] * 100 + 100
            for x in avg_superkmer_size
        ]
        avg_nb_sk_diff_pct = [
            (x - nb_superkmer[0]) / nb_superkmer[0] * 100 + 100 for x in nb_superkmer
        ]

        fig, axes = plt.subplots(4, 2, figsize=(12, 6))

        # Original plots
        axes[0, 0].plot(range1(overhead), overhead, "+--")
        axes[0, 0].plot(range1(overhead), [limit_overhead] * len(overhead))
        axes[0, 0].set_title("base per $k$-mer")
        axes[0, 0].set_xlabel("number of hash function")
        axes[0, 0].set_ylabel("#base in representation / |read|")

        axes[1, 0].plot(range1(avg_superkmer_size), avg_superkmer_size, "+--")
        axes[1, 0].plot(
            range1(avg_superkmer_size),
            [limit_avg_superkmer_size] * len(avg_superkmer_size),
        )
        axes[1, 0].set_title("average superkmer size")
        axes[1, 0].set_xlabel("number of hash function")
        axes[1, 0].set_ylabel("average superkmer size")

        axes[2, 0].plot(range1(nb_superkmer), nb_superkmer, "+--")
        axes[2, 0].plot(
            range1(nb_superkmer),
            [limit_nb_superkmer] * len(nb_superkmer),
        )
        axes[2, 0].set_title("Number of superkmer")
        axes[2, 0].set_xlabel("number of hash function")
        axes[2, 0].set_ylabel("number of superkmer")

        axes[3, 0].plot(range1(density), density, "+--")
        axes[3, 0].set_title("Density")
        axes[3, 0].set_xlabel("number of hash function")
        axes[3, 0].set_ylabel("density")

        # % difference plots
        axes[0, 1].plot(range1(overhead_diff_pct), overhead_diff_pct, "+--")
        axes[0, 1].set_title("Change of base per $k$-mer")
        axes[0, 1].set_xlabel("number of hash function")
        axes[0, 1].set_ylabel("Variation (%)")

        axes[1, 1].plot(range1(avg_size_diff_pct), avg_size_diff_pct, "+--")
        axes[1, 1].set_title("Change of average superkmer size")
        axes[1, 1].set_xlabel("number of hash function")
        axes[1, 1].set_ylabel("Variation (%)")

        axes[2, 1].plot(range1(avg_nb_sk_diff_pct), avg_nb_sk_diff_pct, "+--")
        axes[2, 1].set_title("Change of number of superkmer")
        axes[2, 1].set_xlabel("number of hash function")
        axes[2, 1].set_ylabel("Variation (%)")

    plt.tight_layout()
    plt.show()
