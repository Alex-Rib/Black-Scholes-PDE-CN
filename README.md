# Crank-Nicolson Scheme for the Black-Scholes PDE

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Finance](https://img.shields.io/badge/Finance-Derivatives-green)
![Tests](https://img.shields.io/badge/Tests-pytest-purple)
![Lint](https://img.shields.io/badge/Lint-ruff-orange)
![Status](https://img.shields.io/badge/Status-Educational-orange)

## 📊 Description

Implementation of a **Crank-Nicolson finite-difference scheme** to solve the Black-Scholes partial differential equation (PDE) in log-price and price a European call.

## 🎯 Objectives

- Solve the Black-Scholes PDE by finite differences in log-price coordinates.
- Analyse the error convergence as a function of the time and space steps.
- Measure execution times for different grid resolutions.

## 📐 Mathematical Model

Setting $x = \ln(S)$, the option price $u(t, x)$ satisfies the PDE:

$$\frac{\partial u}{\partial t} + \left(r - \frac{1}{2}\sigma^2\right)\frac{\partial u}{\partial x} + \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2} - r u = 0$$

### Terminal and boundary conditions (Dirichlet)

- **Terminal condition**:
  $$u(T, x) = \max\left(e^x - K, 0\right)$$
- **Lower boundary** ($x \to x_{min}$): $u = 0$
- **Upper boundary** ($x \to x_{max}$): $u = e^x - K e^{-r(T-t)}$

The space grid is centred on $\ln(S_0)$ and its half-width is set from a Gaussian quantile of the terminal log-price distribution, so the grid node $M/2$ corresponds exactly to $S_0$ ($M$ even).

## 🔧 Numerical Method

### Discretisation

- **Time step**: $h = T/N$
- **Space step**: $\delta = (x_{max} - x_{min})/M$

### Crank-Nicolson scheme

The scheme averages the explicit and implicit discretisations (second order in both time and space) and is **unconditionally stable**. At each time step a tridiagonal system is solved:

$$A\,u_{i-1}^{n-1} + B\,u_i^{n-1} + C\,u_{i+1}^{n-1} = \alpha\,u_{i-1}^{n} + \beta\,u_i^{n} + \gamma\,u_{i+1}^{n}$$

with implicit-side coefficients:

$$A = \frac{h(r - 0.5\sigma^2)}{4\delta} - \frac{h\sigma^2}{4\delta^2}, \quad
B = 1 + \frac{rh}{2} + \frac{h\sigma^2}{2\delta^2}, \quad
C = -\left(\frac{h(r - 0.5\sigma^2)}{4\delta} + \frac{h\sigma^2}{4\delta^2}\right)$$

and explicit-side coefficients:

$$\alpha = \frac{h\sigma^2}{4\delta^2} - \frac{h(r - 0.5\sigma^2)}{4\delta}, \quad
\beta = 1 - \frac{rh}{2} - \frac{h\sigma^2}{2\delta^2}, \quad
\gamma = \frac{h(r - 0.5\sigma^2)}{4\delta} + \frac{h\sigma^2}{4\delta^2}$$

The tridiagonal system is solved with SciPy's `solve_banded`; the right-hand side is assembled with vectorised NumPy slicing.

## 📁 Project Structure

```
.
├── black_scholes.py     # Closed-form Black-Scholes call price
├── cn_scheme.py         # OptionParams, GridConfig, CrankNicolsonScheme
├── run_experiments.py   # Convergence study + log-log error plot
└── tests/
    └── Test_cn_scheme.py
```

## 📈 Experiments

`run_experiments.py`:
- computes numerical prices for several $(M, N)$ pairs,
- compares them to the analytical Black-Scholes price,
- plots the error in log-log scale (one curve per $M$).

## 🚀 Usage

```bash
python run_experiments.py
```

Minimal pricing example:

```python
from cn_scheme import CrankNicolsonScheme, GridConfig, OptionParams

params = OptionParams(S0=100, K=100, T=1.0, r=0.05, sigma=0.2)
grid = GridConfig(M=800, N=2000)
price = CrankNicolsonScheme(params, grid).price()
```

## ✅ Tests


Six unit tests cover: convergence to the closed-form price, second-order time accuracy with few time steps, error reduction under grid refinement, unconditional stability (no blow-up at very large time steps), no-arbitrage lower bound, and grid validation.

## 🔗 Related Repositories

Companion repositories implementing the **explicit scheme** (with stability analysis) and the **implicit scheme** (Thomas algorithm vs SciPy `solve_banded`) for the same PDE.

## 👨‍💻 Author

Alexandre R. - Université Paris Cité