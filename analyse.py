import matplotlib.pyplot as plt

with open("data.txt", "r") as fichier:
    contents = fichier.readlines()[0]
    contents = eval(contents)  # TODO dangerous
    print(contents)
    overhead = [x[0] for x in contents]
    avg_size = [x[1] for x in contents]

    # Compute percentage differences from the first element
    overhead_diff_pct = [(x - overhead[0]) / overhead[0] * 100 + 100 for x in overhead]
    avg_size_diff_pct = [(x - avg_size[0]) / avg_size[0] * 100 + 100 for x in avg_size]

    fig, axes = plt.subplots(2, 2, figsize=(12, 6))

    # Original plots
    axes[0, 0].plot(range(1, len(overhead) + 1), overhead)
    axes[0, 0].set_title("base per $k$-mer")
    axes[0, 0].set_xlabel("number of hash function")
    axes[0, 0].set_ylabel("#base in representation / |read|")

    axes[1, 0].plot(range(1, len(avg_size) + 1), avg_size)
    axes[1, 0].set_title("average superkmer size")
    axes[1, 0].set_xlabel("number of hash function")
    axes[1, 0].set_ylabel("average superkmer size")

    # % difference plots
    axes[0, 1].plot(range(len(overhead_diff_pct)), overhead_diff_pct)
    axes[0, 1].set_title("Change of base per $k$-mer")
    axes[0, 1].set_xlabel("number of hash function")
    axes[0, 1].set_ylabel("% Size with one hash function")

    axes[1, 1].plot(range(len(avg_size_diff_pct)), avg_size_diff_pct)
    axes[1, 1].set_title("Change of average superkmer size")
    axes[1, 1].set_xlabel("number of hash function")
    axes[1, 1].set_ylabel("% Size with one hash function")

    plt.tight_layout()
    plt.show()
