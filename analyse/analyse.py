import matplotlib.pyplot as plt
import json

filenames = ["data_canonical.json", "data_non_canonical.json"]


def range1(collection):
    return range(1, len(collection) + 1)


def roundvec(v, precision):
    return [round(x, precision) for x in v]


def plot_single_k(content, outfilename):
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
    density = [x / (size_read - m + 1) for x in nb_superkmer]
    limit_density = 1 / (k - m + 1)

    # Compute percentage differences from the first element
    overhead_diff_pct = [(x - overhead[0]) / overhead[0] * 100 + 100 for x in overhead]
    avg_size_diff_pct = [
        (x - avg_superkmer_size[0]) / avg_superkmer_size[0] * 100 + 100
        for x in avg_superkmer_size
    ]
    avg_nb_sk_diff_pct = [
        (x - nb_superkmer[0]) / nb_superkmer[0] * 100 + 100 for x in nb_superkmer
    ]

    fig, axes = plt.subplots(4, 2, figsize=(12, 6))
    # fig.canvas.manager.set_window_title(f"{filename=}, {k=}, {m=}")

    # let's round everything
    precision = 3
    limit_overhead = round(limit_overhead, precision)
    overhead = roundvec(overhead, precision)
    limit_avg_superkmer_size = round(limit_avg_superkmer_size, precision)
    avg_superkmer_size = roundvec(avg_superkmer_size, precision)
    limit_nb_superkmer = round(limit_nb_superkmer, precision)
    nb_superkmer = roundvec(nb_superkmer, precision)

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
    # axes[2, 0].set_ylim(bottom=0)

    axes[3, 0].plot(range1(density), density, "+--")
    axes[3, 0].plot(range1(density), [limit_density] * len(density))
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
    # plt.show()
    plt.savefig(outfilename)
    plt.clf()


def plot_vs_w(content_all_k, outfilename):
    ws = []
    densities = [[] for _ in range(16)]
    for k in content_all_k:
        content = content_all_k[k]
        k = int(k)
        size_read = content["size_read"]
        k = content["k"]
        m = content["m"]
        nb_superkmer = content["nb_superkmer"]

        # compute the density
        for i in range(16):
            densities[i].append(nb_superkmer[i] / (size_read - m + 1))

        ws.append(k - m + 1)
        # densities.append(density)
    # print(densities)
    for i in range(len(densities)):
        densities[i] = [x for _, x in sorted(zip(ws, densities[i]))]
    ws.sort()
    for i in range(16):
        plt.plot(ws, densities[i], "+--", label=f"{i + 1} hash functions")
    plt.legend()
    plt.xlabel("w")
    plt.ylabel("density")

    plt.tight_layout()
    plt.savefig(outfilename)
    plt.clf()


def main():
    # plot density figure
    for filename in filenames:
        with open(filename, "r") as fichier:
            content_all_k = json.load(fichier)
            content_all_k = content_all_k["data"]
            outfilename = "analyse/plot_density" + filename.removesuffix(".json")
            plot_vs_w(content_all_k, outfilename)

    # plot single k figures
    for filename in filenames:
        with open(filename, "r") as fichier:
            content_all_k = json.load(fichier)
            content_all_k = content_all_k["data"]
            for k in content_all_k:
                content = content_all_k[k]
                # print(content)
                k = int(k)
                outfilename = (
                    "analyse/analyse_" + filename.removesuffix(".json") + f"_{k=}.png"
                )
                plot_single_k(content, outfilename)


if __name__ == "__main__":
    main()
