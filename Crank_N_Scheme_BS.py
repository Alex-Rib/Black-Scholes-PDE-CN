####################################
# import des librairies necessaires
####################################
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.linalg import solve_banded


#####################
# paramètres globaux
#####################

# paramètres de l'option
S0 = 100  
K = 100   
T = 1.0   
r = 0.05  
sigma = 0.2  


# paramètres pour le schéma implicite 
M = 800
N = 20000
z = norm.ppf(0.99)
l = np.abs((r - 0.5 * sigma**2)) * T + (z * sigma * np.sqrt(T)) 
delta_min = np.log(S0) - l
delta_max = np.log(S0) + l
delta = (delta_max - delta_min) / M 
h = T / N  
delta_values = np.linspace(delta_min, delta_max, M+1) 

###########################
# définition des fonctions
###########################

def Crank_Nicolson_BS_log(S0, K, T, r, sigma, M, N,h,delta,delta_values,delta_max):
    S_values = np.exp(delta_values)
    u_price = np.maximum(S_values - K, 0.0)

    A = ((h*(r - 0.5*sigma**2))/(4*delta) - (h*sigma**2)/(4*delta**2)) # coefficient de u_{i,j-1}
    B = ( 1 + (0.5*(r*h)) + (h*sigma**2)/(2*delta**2))                     # coefficient de u_{i,j}
    C = -((h*(r-0.5*sigma**2))/(4*delta) + (h*sigma**2)/(4*delta**2)) # coefficient de u_{i,j+1}

    ALPHA = ((h*sigma**2)/(4*delta**2)-(h*(r - 0.5*sigma**2))/(4*delta)) # coefficient de u_{i+1,j-1}
    BETA = ( 1 - (0.5*(r*h)) - (h*sigma**2)/(2*delta**2))                     # coefficient de u_{i+1,j}
    GAMMA = ((h*(r-0.5*sigma**2))/(4*delta) + (h*sigma**2)/(4*delta**2)) # coefficient de u_{i+1,j+1}

    MATRICE_G = np.zeros((3,M-1))
    MATRICE_G[2,:-1] = A  # sous-diagonale
    MATRICE_G[1,:] = B    # diagonale
    MATRICE_G[0,1:] = C   # sur-diagonale

    u_droite = np.zeros(M-1)
    exp_delta = np.exp(delta_max)

    for i in range(N,0,-1):
        t_old = i*h
        t_new = (i-1) * h

        u_low_old = u_low_new = 0.0  
        u_high_old = exp_delta - K * np.exp(-r * (T - t_old))  
        u_high_new = exp_delta - K * np.exp(-r * (T - t_new))

        u_droite[:] = u_price[1:M] * BETA
        u_droite[1:] += u_price[1: M-1] * ALPHA
        u_droite[:-1] += u_price[2:M] * GAMMA

        # conditions aux limites (Dirichlet)
        u_droite[0] += ALPHA * u_low_old - A * u_low_new
        u_droite[-1] += GAMMA * u_high_old - C * u_high_new

        u_mid = solve_banded((1,1), MATRICE_G, u_droite)
        u_price[0] = u_low_new
        u_price[1:M] = u_mid
        u_price[M] =  u_high_new

    return u_price[M//2]  # on retourne le prix à S0 (au milieu de la grille)


def call_bs(S, K, T, r, sigma):
    d1 = (np.log(S/K) + (r + 0.5 * sigma**2)*T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    price = S * norm.cdf(d1) - K * np.exp(-r*T)* norm.cdf(d2)
    return price



def main():
#############################################################################
# Analyse de la convergence du schema de Crank-Nicolson en fonction de M et N
#############################################################################

    M_values = [100,400,800,1600,3200,6400,12800]
    N_values = [100,200,400,600,800,1000,2000,4000]
    Matrice_prix = np.zeros((len(M_values), len(N_values)))
    Matrice_temps = np.zeros((len(M_values), len(N_values)))
    Matrice_erreur = np.zeros((len(M_values), len(N_values)))
    bs_price = call_bs(S0,K,T,r,sigma)
    data = []
    for i, m in enumerate(M_values):
        delta_m = (delta_max - delta_min)/m 
        delta_values_m = np.linspace(delta_min, delta_max, m+1)
        for j, n in enumerate(N_values):
            h_n = T/n 
            start_time = time.time()
            Matrice_prix[i,j] = Crank_Nicolson_BS_log(S0, K, T, r, sigma, m, n, h_n, delta_m, delta_values_m, delta_max)
            end_time = time.time()
            Matrice_temps[i,j] = end_time - start_time
            Matrice_erreur[i,j] = np.abs(Matrice_prix[i,j] - bs_price)
            data.append({'M': m,'N': n,'Prix': Matrice_prix[i,j], 'Prix BS': bs_price,'Erreur': Matrice_erreur[i,j], 'Temps': Matrice_temps[i,j]})
    data_frame = pd.DataFrame(data)
    print(data_frame)

    plt.figure(figsize=(12,6))
    for i, m in enumerate(M_values):
        plt.plot(N_values, Matrice_erreur[i,:], marker='o', label=f'M={m}')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Nombre de points temporels N")
    plt.ylabel('Erreur')
    plt.title("Convergence de l'erreur du schéma de Crank- Nicolson")
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()


if __name__ == "__main__":
    main()

        


