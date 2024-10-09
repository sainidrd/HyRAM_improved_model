# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 12:27:40 2022

@author: saini
"""


import numpy as np

Tcrit = np.array([33.2, 126.2, 154.581])  # Critical temperature in Kelvin
Pcrit = np.array([1.297e6, 3.394e6, 5.043e6])  # Critical pressure in Pa
omega = np.array([-0.219, 0.039, 0.0222])  # Acentric factor
R = 8.3145  # Gas constant, J/(mol*K)
MW_H2 = 2.01588e-3
MW_N2 = 0.02896546
Ka_New = [[0.0, 0.0436, 0.433], [0.0436, 0.0, -0.0821], [0.433, -0.0821, 0.0]]
Kb_New = [[0.0, -0.1611, 0.0382], [-0.1611, 0.0, -0.2633], [0.0382, -0.2633, 0.0]]



# Peng-Robinson coefficients calculation for H2, N2, O2 mixture
def PR_H2_N2_Coeff(T, P, x_H2):
    x_Air = 1 - x_H2
    x_N2 = (0.767)*x_Air
    x_O2 = (0.233)*x_Air
    X = np.array([x_H2, x_N2, x_O2])  # Composition of H2, N2, O2

    a = np.zeros(3)
    b = np.zeros(3)
    Tr = np.zeros(3)
    m = np.zeros(3)
    alpha = np.zeros(3)

    for i in range(3):
        Tr[i] = T / Tcrit[i]
        m[i] = 0.37464 + 1.54226 * omega[i] - 0.26992 * omega[i] ** 2
        alpha[i] = (1 + m[i] * (1 - np.sqrt(Tr[i]))) ** 2
        a[i] = 0.45724 * alpha[i] * (R * Tcrit[i]) ** 2 / Pcrit[i]
        b[i] = 0.07780 * (R * Tcrit[i]) / Pcrit[i]

    # Mixture coefficients
    a_m = np.zeros((3, 3))
    b_m = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            a_m[i][j] = (1 - Ka_New[i][j]) * (X[i] * X[j]) * np.sqrt(a[i] * a[j])
            b_m[i][j] = (1 - Kb_New[i][j]) * (X[i] * X[j]) * ((b[i] + b[j]) / 2)

    A = np.sum(a_m) * P / (R * T) ** 2
    B = np.sum(b_m) * P / (R * T)

    return A, B, alpha, Tr, m


def PR_H2_N2_Coeff_ab(T, P, x_H2):
    x_Air = 1 - x_H2
    X = np.array([x_H2, x_Air, 0.0])  # Composition of H2, N2, O2
    x_N2 = (0.767)*x_Air
    x_O2 = (0.233)*x_Air
    X = np.array([x_H2, x_N2, x_O2])  # Composition of H2, N2, O2
    a = np.zeros(3)
    b = np.zeros(3)
    Tr = np.zeros(3)
    m = np.zeros(3)
    alpha = np.zeros(3)

    for i in range(3):
        Tr[i] = T / Tcrit[i]
        m[i] = 0.37464 + 1.54226 * omega[i] - 0.26992 * omega[i] ** 2
        alpha[i] = (1 + m[i] * (1 - np.sqrt(Tr[i]))) ** 2
        a[i] = 0.45724 * alpha[i] * (R * Tcrit[i]) ** 2 / Pcrit[i]
        b[i] = 0.07780 * (R * Tcrit[i]) / Pcrit[i]
        
    # Mixture coefficients
    a_m = np.zeros((3, 3))
    b_m = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            a_m[i][j] = (1 - Ka_New[i][j]) * (X[i] * X[j]) * np.sqrt(a[i] * a[j])
            b_m[i][j] = (1 - Kb_New[i][j]) * (X[i] * X[j]) * ((b[i] + b[j]) / 2)

    return np.sum(a_m), np.sum(b_m)


def PR_H2_N2_Coeff_ab_ind(T, P, x_H2):
    x_Air = 1 - x_H2
    X = np.array([x_H2, x_Air, 0.0])  # Composition of H2, N2, O2
    x_N2 = (0.767)*x_Air
    x_O2 = (0.233)*x_Air
    X = np.array([x_H2, x_N2, x_O2])  # Composition of H2, N2, O2
    a = np.zeros(3)
    b = np.zeros(3)
    Tr = np.zeros(3)
    m = np.zeros(3)
    alpha = np.zeros(3)

    for i in range(3):
        Tr[i] = T / Tcrit[i]
        m[i] = 0.37464 + 1.54226 * omega[i] - 0.26992 * omega[i] ** 2
        alpha[i] = (1 + m[i] * (1 - np.sqrt(Tr[i]))) ** 2
        a[i] = 0.45724 * alpha[i] * (R * Tcrit[i]) ** 2 / Pcrit[i]
        b[i] = 0.07780 * (R * Tcrit[i]) / Pcrit[i]
    return a, b

def find_real_roots(r):
    real_roots = [r[i].real for i in range(len(r)) if r[i].imag == 0.0 and r[i].real > 0.0]
    return real_roots


# Fugacity coefficient calculation
def fugacity(Z, A, B):
    N = len(Z)
    phi = np.zeros(N)
    for i in range(N):
        sq = np.sqrt(2)
        phi[i] = np.exp(Z[i]-1-np.log(Z[i]-B)-A/B/2/sq*np.log((Z[i] + (1+sq)*B)/(Z[i] + (1-sq)*B)))
    return phi


# Sorting and selecting appropriate root
def Sort_roots(real_roots, A, B):
    if len(real_roots) > 1:
        real_roots = [real_roots[0]]  # Select the first root
    phi = fugacity(real_roots, A, B)
    return real_roots, phi

# Main Peng-Robinson function
def PR_H2_N2_Main(T, P, x_H2):
    A, B, _, _, _ = PR_H2_N2_Coeff(T, P, x_H2)

    # Polynomial coefficients for Peng-Robinson EOS
    Z_poly = [1, -(1 - B), A - 3 * B ** 2 - 2 * B, -(A * B - B ** 2 - B ** 3)]
    roots = np.roots(Z_poly)
    real_roots = find_real_roots(roots)
    Z, phi_Z = Sort_roots(real_roots, A, B)
    return Z, phi_Z


###################################################################################################
# Density calculation from pressure and temperature
def find_rho_PT_Mixture(p, T, x_H2):
    mol_mass = (x_H2 * MW_H2 + (1 - x_H2) * MW_N2)
    Z, _ = PR_H2_N2_Main(T, p, x_H2)
    try:
        rho_val = (p/(R*T))*(mol_mass/Z) 
    except:
        rho_val = (p/(R*T))*(mol_mass/Z[0]) 
    return rho_val

# Temperature calculation from density and pressure
def find_T_rhoP_Mixture(rho, p, xH2):
    Y_ = xH2 * MW_H2 / (xH2 * MW_H2 + (1 - xH2) * MW_N2)
    mol_mass = MW_N2 * MW_H2 / (Y_ * (MW_N2 - MW_H2) + MW_H2)

    T = (p * mol_mass) / (rho * R)
    T1 = T - 0.1 * T #10% tolerance
    T2 = T + 0.1 * T

    T_final = regulaFalsi_T(T1, T2, p, rho, xH2)
    return T_final


# Derivative of temperature with respect to density
def derivative_T_rhoP_Mixture(rho, p, xH2):
    T = find_T_rhoP_Mixture(rho, p, xH2)
    a, b = PR_H2_N2_Coeff_ab(T, p, xH2)
    Y_ = xH2 * MW_H2 / (xH2 * MW_H2 + (1 - xH2) * MW_N2)
    mol_mass = MW_N2 * MW_H2 / (Y_ * (MW_N2 - MW_H2) + MW_H2)

    dT_drho = 2 * a * (mol_mass / R) * ((mol_mass ** 2 - (b * rho) ** 2) / (mol_mass ** 2 + 2 * mol_mass * b * rho - (b * rho) ** 2) ** 2) - \
              (mol_mass * T / rho) * (1 / (mol_mass - b * rho))

    return dT_drho


def derivative_T_Y_Mixture(rho, p, xH2):
    T = find_T_rhoP_Mixture(rho, p, xH2)
    a, b = PR_H2_N2_Coeff_ab_ind(T, p, xH2)

    Y_ = xH2 * MW_H2 / (xH2 * MW_H2 + (1 - xH2) * MW_N2)
    mol = MW_N2 * MW_H2 / (Y_ * (MW_N2 - MW_H2) + MW_H2)

    dmoldY = ((MW_H2 - MW_N2) * MW_H2 * MW_N2) / (MW_H2 + Y_ * (MW_N2 - MW_H2)) ** 2
   # PR coefficients for the mixture
    am, bm = PR_H2_N2_Coeff_ab(T, p, xH2)
    a, b = PR_H2_N2_Coeff_ab_ind(T, p, xH2)
    F_N2 = 0.767
    F_O2 = 0.237

   # Calculate dadY and dbdY terms
    dadY = (
       2 * xH2 * a[0] * (1 - Ka_New[0][0])
       - 2 * F_N2 * (1 - xH2) * a[1] * (1 - Ka_New[1][1])
       - 2 * a[2] * F_O2**2 * (1 - xH2) * (1 - Ka_New[2][2])
       + 2 * (1 - Ka_New[0][1]) * F_N2 * (1 - 2 * xH2) * (a[0] * a[1]) ** 0.5
       + 2 * (1 - Ka_New[0][2]) * F_O2 * (1 - 2 * xH2) * (a[0] * a[2]) ** 0.5
       - 4 * (1 - Ka_New[1][2]) * F_N2 * F_O2 * (1 - xH2) * (a[1] * a[2]) ** 0.5
       ) * (MW_H2 * MW_N2) / (MW_H2 + Y_ * (MW_N2 - MW_H2)) ** 2

    dbdY = (
       2 * xH2 * b[0] * (1 - Kb_New[0][0])
       - 2 * F_N2 * (1 - xH2) * b[1] * (1 - Kb_New[1][1])
       - 2 * F_O2 * (1 - xH2) * (1 - Kb_New[2][2]) * b[2]
       + (1 - Kb_New[0][1]) * F_N2 * (1 - 2 * xH2) * (b[0] + b[1])
       + (1 - Kb_New[0][2]) * F_O2 * (1 - 2 * xH2) * (b[0] + b[2])
       - (1 - Kb_New[1][2]) * F_N2 * F_O2 * (1 - xH2) * (b[1] + b[2])
       ) * (MW_H2 * MW_N2) / (MW_H2 + Y_ * (MW_N2 - MW_H2)) ** 2

   # Final derivative calculation
    dTdY = (
       (T / (mol - bm * rho)) * (dmoldY - rho * dbdY)
       + (rho / R) * ((mol - bm * rho) / (mol**2 + 2 * mol * bm * rho - (bm * rho)**2))
       * (dadY - (2 * am / (mol**2 + 2 * mol * bm * rho - (bm * rho)**2))
       * ((mol + bm * rho) * dmoldY + (mol - bm * rho) * rho * dbdY))
       )

    return dTdY

##############################################################################################################################
def func_TPrho(T, p, rho, xH2):
    # Calculate molar mass of the mixture
    Y_ = (xH2 * MW_H2) / (xH2 * MW_H2 + (1 - xH2) * MW_N2)
    mol_mass = (MW_N2 * MW_H2) / (Y_ * (MW_N2 - MW_H2) + MW_H2)

    # Get compressibility factor (Z) from PR equation
    Z, _ = PR_H2_N2_Main(T, p, xH2)

    # Calculate X based on Z
    try:
        X = p / Z - rho * R * T / mol_mass
    except TypeError:
        X = p / Z[0] - rho * R * T / mol_mass

    return X



  
def regulaFalsi_T( a , b, p, rho, xH2):
    MAX_ITER = 1000000
    if func_TPrho(a, p, rho, xH2) * func_TPrho(b, p, rho, xH2) >= 0:
        print("You have not assumed right a and b")
        return -1
     
    c = a # Initialize result
     
    for i in range(MAX_ITER):
        # Find the point that touches x axis
        c = (a * func_TPrho(b, p, rho, xH2) - b * func_TPrho(a, p, rho, xH2))/ (func_TPrho(b, p, rho, xH2) - func_TPrho(a, p, rho, xH2))
         
        # Check if the above found point is root
        if func_TPrho(c, p, rho, xH2) == 0:
            break
         
        # Decide the side to repeat the steps
        elif func_TPrho(c, p, rho, xH2) * func_TPrho(a, p, rho, xH2) < 0:
            b = c
        else:
            a = c
        return c
