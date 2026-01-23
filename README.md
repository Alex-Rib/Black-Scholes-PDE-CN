# Schéma de Crank–Nicolson pour l'EDP de Black–Scholes

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![Finance](https://img.shields.io/badge/Finance-Derivatives-green)
![Status](https://img.shields.io/badge/Status-Educational-orange)

## 📊 Description

Implémentation d'un **schéma de Crank–Nicolson** pour résoudre l’EDP de Black–Scholes (call européen) dans la variable $x = \ln(S)$. Le script étudie la **convergence** du schéma en fonction de $M$ (espace) et $N$ (temps), et compare les prix numériques au prix fermé de Black–Scholes.

## 🎯 Objectifs

- Résoudre l’EDP de Black–Scholes par différences finies en log‑prix.
- Analyser la convergence de l’erreur en fonction des pas temporel et spatial.
- Mesurer les temps d’exécution pour différentes grilles.

## 📐 Modèle Mathématique

Après le changement de variable $x = \ln(S)$, l’EDP s’écrit :

$$\frac{\partial u}{\partial t} + \left(r - \frac{1}{2}\sigma^2\right)\frac{\partial u}{\partial x} + \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2} - r u = 0$$

### Conditions aux bornes (Dirichlet) et terminale

- **Condition terminale** :
  $$u(T, x) = \max(e^x - K, 0)$$
- **Borne inférieure** ($x \to x_{min}$) : $u = 0$
- **Borne supérieure** ($x \to x_{max}$) : $u = e^x - K e^{-r(T-t)}$

## 🔧 Méthode Numérique

### Discrétisation

- **Pas de temps** : $h = T/N$
- **Pas d’espace** : $\delta = (x_{max} - x_{min})/M$

### Schéma de Crank–Nicolson

Le schéma combine l’explicite et l’implicite (ordre 2 en temps) :

$$A\,u_{i-1}^{n-1} + B\,u_i^{n-1} + C\,u_{i+1}^{n-1} = \alpha\,u_{i-1}^{n} + \beta\,u_i^{n} + \gamma\,u_{i+1}^{n}$$

avec les coefficients utilisés dans le script :

$$A = \frac{h(r - 0.5\sigma^2)}{4\delta} - \frac{h\sigma^2}{4\delta^2}, \quad
B = 1 + \frac{rh}{2} + \frac{h\sigma^2}{2\delta^2}, \quad
C = -\left(\frac{h(r - 0.5\sigma^2)}{4\delta} + \frac{h\sigma^2}{4\delta^2}\right)$$

$$\alpha = \frac{h\sigma^2}{4\delta^2} - \frac{h(r - 0.5\sigma^2)}{4\delta}, \quad
\beta = 1 - \frac{rh}{2} - \frac{h\sigma^2}{2\delta^2}, \quad
\gamma = \frac{h(r - 0.5\sigma^2)}{4\delta} + \frac{h\sigma^2}{4\delta^2}$$

La résolution du système tridiagonal est faite via `solve_banded` (SciPy).

## 📈 Expériences

Le script :
- calcule les prix numériques pour plusieurs couples $(M, N)$ ;
- compare au prix analytique de Black–Scholes ;
- trace l’erreur en échelle log–log.

## 🚀 Utilisation

```bash
python Crank_N_Scheme_BS.py
```

## 📦 Dépendances

```bash
pip install numpy pandas matplotlib scipy
```

## 👨‍💻 Auteur

Alexandre R. - Master ISIFAR, Université Paris Cité
