"""Crank-Nicolson finite-difference scheme for the Black-Scholes PDE in log-price."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_banded
from scipy.stats import norm


@dataclass(frozen=True)
class OptionParams:
    """European call option parameters."""

    S0: float = 100.0
    K: float = 100.0
    T: float = 1.0
    r: float = 0.05
    sigma: float = 0.2


@dataclass(frozen=True)
class GridConfig:
    """Space-time grid for the Crank-Nicolson scheme.

    The space grid is in log-price, centred on log(S0) and sized from a
    Gaussian quantile of the terminal log-price distribution, so that the
    index M // 2 corresponds exactly to S0 (M must be even).
    """

    M: int = 800
    N: int = 20_000
    confidence: float = 0.99

    def __post_init__(self) -> None:
        if self.M % 2 != 0:
            raise ValueError("M must be even so that the grid is centred on S0.")

    def space_grid(self, params: OptionParams) -> np.ndarray:
        """Return the log-price grid (M + 1 points centred on log(S0))."""
        z = norm.ppf(self.confidence)
        half_width = abs(params.r - 0.5 * params.sigma**2) * params.T + z * params.sigma * np.sqrt(
            params.T
        )
        x_min = np.log(params.S0) - half_width
        x_max = np.log(params.S0) + half_width
        return np.linspace(x_min, x_max, self.M + 1)


class CrankNicolsonScheme:
    """Backward-in-time Crank-Nicolson scheme (second order in time and space).

    Each time step solves the tridiagonal system
    A u_{i-1}^{n-1} + B u_i^{n-1} + C u_{i+1}^{n-1}
        = alpha u_{i-1}^n + beta u_i^n + gamma u_{i+1}^n,
    averaging the explicit and implicit discretisations, with Dirichlet
    boundary conditions. The scheme is unconditionally stable.
    """

    def __init__(self, params: OptionParams, grid: GridConfig) -> None:
        self.params = params
        self.grid = grid

    def price(self) -> float:
        """Return the call price at (t = 0, S = S0)."""
        p, g = self.params, self.grid
        x = g.space_grid(p)
        dx = x[1] - x[0]
        h = p.T / g.N
        m = g.M

        u = np.maximum(np.exp(x) - p.K, 0.0)  # terminal payoff

        drift = p.r - 0.5 * p.sigma**2

        # implicit side (left-hand side, time n-1)
        a = h * drift / (4 * dx) - h * p.sigma**2 / (4 * dx**2)
        b = 1 + 0.5 * p.r * h + h * p.sigma**2 / (2 * dx**2)
        c = -(h * drift / (4 * dx) + h * p.sigma**2 / (4 * dx**2))

        # explicit side (right-hand side, time n)
        alpha = h * p.sigma**2 / (4 * dx**2) - h * drift / (4 * dx)
        beta = 1 - 0.5 * p.r * h - h * p.sigma**2 / (2 * dx**2)
        gamma = h * drift / (4 * dx) + h * p.sigma**2 / (4 * dx**2)

        band = np.zeros((3, m - 1))
        band[0, 1:] = c
        band[1, :] = b
        band[2, :-1] = a

        rhs = np.zeros(m - 1)
        x_max_exp = np.exp(x[-1])

        for i in range(g.N, 0, -1):
            t_old = i * h
            t_new = (i - 1) * h

            # Dirichlet boundary conditions at both time levels
            u_low_old = u_low_new = 0.0
            u_up_old = x_max_exp - p.K * np.exp(-p.r * (p.T - t_old))
            u_up_new = x_max_exp - p.K * np.exp(-p.r * (p.T - t_new))

            # explicit combination at time n
            rhs[:] = beta * u[1:m]
            rhs[1:] += alpha * u[1 : m - 1]
            rhs[:-1] += gamma * u[2:m]

            # boundary contributions
            rhs[0] += alpha * u_low_old - a * u_low_new
            rhs[-1] += gamma * u_up_old - c * u_up_new

            u_mid = solve_banded((1, 1), band, rhs)

            u[0] = u_low_new
            u[1:m] = u_mid
            u[m] = u_up_new

        return float(u[m // 2])
