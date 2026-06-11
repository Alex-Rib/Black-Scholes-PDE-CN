"""Unit tests for the Crank-Nicolson Black-Scholes PDE scheme."""

import numpy as np
import pytest

from black_scholes import call_price
from cn_scheme import CrankNicolsonScheme, GridConfig, OptionParams

PARAMS = OptionParams()
BS_PRICE = call_price(PARAMS.S0, PARAMS.K, PARAMS.T, PARAMS.r, PARAMS.sigma)


def test_converges_to_closed_form():
    """CN price converges to the closed-form Black-Scholes price."""
    cn_price = CrankNicolsonScheme(PARAMS, GridConfig(M=800, N=2000)).price()
    assert cn_price == pytest.approx(BS_PRICE, abs=5e-3)


def test_accurate_with_few_time_steps():
    """Second-order accuracy in time: good precision with a coarse time grid."""
    cn_price = CrankNicolsonScheme(PARAMS, GridConfig(M=800, N=200)).price()
    assert cn_price == pytest.approx(BS_PRICE, abs=1e-2)


def test_error_decreases_with_refinement():
    """Refining the grid reduces the error against the closed-form price."""
    coarse = CrankNicolsonScheme(PARAMS, GridConfig(M=100, N=100)).price()
    fine = CrankNicolsonScheme(PARAMS, GridConfig(M=800, N=800)).price()
    assert abs(fine - BS_PRICE) < abs(coarse - BS_PRICE)


def test_stable_with_large_time_steps():
    """Unconditional stability: no blow-up even with very few time steps."""
    price = CrankNicolsonScheme(PARAMS, GridConfig(M=800, N=20)).price()
    assert np.isfinite(price)
    assert 0.0 < price < PARAMS.S0


def test_price_positive_and_above_intrinsic():
    """Call price is positive and above the no-arbitrage lower bound."""
    price = CrankNicolsonScheme(PARAMS, GridConfig(M=400, N=400)).price()
    lower_bound = PARAMS.S0 - PARAMS.K * np.exp(-PARAMS.r * PARAMS.T)
    assert price > 0.0
    assert price >= lower_bound


def test_odd_m_raises():
    """An odd number of space points is rejected (grid not centred on S0)."""
    with pytest.raises(ValueError):
        GridConfig(M=201, N=100)
