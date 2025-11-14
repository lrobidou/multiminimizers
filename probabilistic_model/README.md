These scripts allow to reproduce some simulations and figures from the article.

In particular,
- `random_minimizers_density.py` allows to operate Monte-Carlo simulations of random minimizers under our probabilistic model, and is interested especially in the gap between consecutive selected positions, thus highlighting the connection between density and gaps, as per (one of) our main result of the paper.
- `random_minimizers_empirical_densities_comparison.py` generates long random sequences (still under our model) and compute the standard and the deduplicated densities every 1000 windows, to hightlight how theses quantites evolve over the course of parsing a sequence.
- finally, `random_minimizers_deduplicated_density_theory.py` compute the theoretical value of the deduplicated density for small values of k in the context of random minimizers. It does not evolve Monte-Carlo simulation.

All these scripts are configured with the parameters used in the paper (they can be modified in the first few lines of the scripts) and can be executed with the command `python script.py`.

Dependencies
------------
`scipy`, `Decimal`, `matplotlib`.
