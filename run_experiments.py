"""Convergence experiments for the Crank-Nicolson scheme.

Run with: python run_experiments.py
"""

import time

import matplotlib.pyplot as plt
import numpy as np

from black_scholes import call_price
from cn_scheme import CrankNicolsonScheme, GridConfig, OptionParams


def convergence_study(
    params: OptionParams, m_values: list[int], n_values: list[int]
) -> tuple[np.ndarray, float]:
    """Compute prices and errors over a grid of (M, N) pairs and print a table."""
    errors = np.zeros((len(m_values), len(n_values)))
    bs_price = call_price(params.S0, params.K, params.T, params.r, params.sigma)

    header = f"{'M':>7} {'N':>6} {'Time (s)':>10} {'CN price':>12} {'BS price':>10} {'Error':>12}"
    print(header)
    print("-" * len(header))

    for i, m in enumerate(m_values):
        for j, n in enumerate(n_values):
            grid = GridConfig(M=m, N=n)
            start = time.perf_counter()
            cn_price = CrankNicolsonScheme(params, grid).price()
            elapsed = time.perf_counter() - start
            errors[i, j] = abs(cn_price - bs_price)
            print(
                f"{m:>7} {n:>6} {elapsed:>10.3f} {cn_price:>12.6f} "
                f"{bs_price:>10.6f} {errors[i, j]:>12.2e}"
            )

    return errors, bs_price


def plot_error(m_values: list[int], n_values: list[int], errors: np.ndarray) -> None:
    """Plot the absolute error against N in log-log scale, one curve per M."""
    plt.figure(figsize=(12, 6))
    for i, m in enumerate(m_values):
        plt.plot(n_values, errors[i, :], marker="o", label=f"M = {m}")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Number of time steps N")
    plt.ylabel("Absolute error")
    plt.title("Error convergence of the Crank-Nicolson scheme")
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()


def main() -> None:
    """Run the convergence study and plot the error curves."""
    params = OptionParams()
    m_values = [100, 400, 800, 1600, 3200, 6400, 12800]
    n_values = [100, 200, 400, 600, 800, 1000, 2000, 4000]

    errors, _ = convergence_study(params, m_values, n_values)
    plot_error(m_values, n_values, errors)


if __name__ == "__main__":
    main()
