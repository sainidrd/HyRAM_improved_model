# -*- coding: utf-8 -*-
"""
Created on Tue May 24 17:29:21 2022

@author: saini
"""
import numpy as np
from scipy import optimize
import math as mt
global Tcrit, Pcrit
global Tref, Pref, tol
Tref =  20.3689#33.145
Pref = 101325.0#1296400.0
tol = 20.0

def PR_H2_Coeff(T, P):
    Tcrit = 33.145 #in kelvin
    Pcrit = 1296400.0 # in Pa
    Tr    = T/Tcrit
    omega = -0.219
    R     = 8.3145# 2.016 # kg/kmol Gas Constant
    m     = 0.37464 + 1.54226*omega - 0.26992*omega**2
    alpha = (1+m*(1-np.sqrt(Tr)))**2
    a     = 0.45724*alpha*((R*Tcrit)**2)/Pcrit
    b     = 0.07780*(R*Tcrit)/Pcrit
    A     = (a*P)/(R*T)**2
    B     = (b*P)/(R*T)
    return A,B,alpha,Tr,m

def find_real_roots(r):
    N = len(r)
    #print ("r_new=", r[2].imag)
    for i in range(N):
        rr=[]
        if r[i].real>0:
            if r[i].imag==0:
               RR = r[i]
            #print ("r_new=", r_new)
    rr.append(RR.real)
    #print ("rr=",rr)
    return rr

def fugacity(Z, A, B):
    N = len(Z)
    phi = np.zeros(N)
    for i in range(N):
        sq = np.sqrt(2)
        phi[i] = np.exp(Z[i]-1-np.log(Z[i]-B)-A/B/2/sq*np.log((Z[i] + (1+sq)*B)/(Z[i] + (1-sq)*B)))
    return phi

def Sort_roots(real_roots, A, B):
    if len(real_roots) > 1:
        print ("More than one roots exist for PR equations")
        XX  = np.where(real_roots > B )[0] # index where roots are > B
        phi = fugacity(real_roots[XX], A, B) # find fugacity of those roots
        if len(phi) > 1:
            XXX = np.where(phi == min(phi))[0]
            real_roots = real_roots[XXX]
            phi = phi[XXX]
    elif len(real_roots)==1:
        print ("Just one root exist for PR equation")
        phi = fugacity(real_roots, A, B)
    return real_roots, phi
    

def PR_H2_Main(T, P):
    
    #Calculate Coeffs for Peng-Robinson Polynomial
    A,B,alpha,Tr,m= PR_H2_Coeff(T,P)
    
    # Setup the Polynomial Coeff
    a1 = 1.0
    a2 = -(1.-B)
    a3 = (A-3.*B*B-2.*B)
    a4 = -(A*B-B*B-B*B*B)
    
    #Determine the roots of cubic polynomial
    Z_poly = [a1, a2, a3, a4]
    roots  = np.roots(Z_poly)
    #print (roots[0])
    #print (roots[1])
    #print (roots[2])
    real_roots = find_real_roots(roots)
    ZZ , phi_ZZ = Sort_roots(real_roots, A, B)
    return ZZ, phi_ZZ

def Cal_Therm(P,T):
    cp = 11909.0 #10587.502897752502 #13.12*10**3 #J/kg K
    R = 8.3145 #2.016
    # Ideal gas contribution
    #H_ideal = cp*(T-Tref)
    H_ideal = (-0.0018/4)*(T**4-Tref**4) + (0.6787/3)*(T**3-Tref**3)\
              + (-57.322/2)*(T**2-Tref**2)+ 11909.0*(T-Tref)#cp*(T-Tref)
    #S_ideal = cp*np.log(T/Tref) - R*np.log(P/Pref)
    S_ideal = (-0.0018/3)*(T**3-Tref**3)+ (0.6787/2)*(T**2-Tref**2) \
             + (-57.322)*(T-Tref) + 11909.0*np.log(T/Tref) - R*np.log(P/Pref) 
    
    # Departure from ideal gas behaviour
    Z, phi = PR_H2_Main(T, P)
    #print ("Z in main=", Z)
    A, B, alpha, Tr, m = PR_H2_Coeff(T, P)
    H_dep = Cal_H_dep(Z, A, B, T, alpha, Tr, R, m)
    S_dep = Cal_S_dep(Z, A, B, T, alpha, Tr, R, m)
    Z_ref, phi_ref = PR_H2_Main(Tref, Pref)
    #print ("Z_ref in main=", Z_ref)
    H_ref = Cal_H_dep(Z_ref, A, B, Tref, alpha, Tr, R, m)
    S_ref = Cal_S_dep(Z_ref, A, B, Tref, alpha, Tr, R, m)
    #S_ref = H_ref/Tref - R*np.log(phi_ref)
    #S_dep = H_dep/T - R*np.log(phi)
    H = H_dep + H_ideal - H_ref
    S = S_dep + S_ideal - S_ref 
    return H, H_ref, S, S_ref


def Cal_H_dep(Z, A, B, T, alpha, Tr, R, m):
    #print ("Z=",Z[0])
    sqrt2 = np.sqrt(2)
    #print ("log function=", Z[0]+(1-sqrt2)*B)
    #print ("B=", B)
    #print ("Z[0]=", Z[0])
    Hdep = R*T*(Z[0]-1.-A/B/2/sqrt2*np.log((Z[0]+(1+sqrt2)*B)/(Z[0]+(1-sqrt2)*B))*(1+m*np.sqrt(Tr)/np.sqrt(alpha)))
    return Hdep
    
def Cal_S_dep(Z, A, B, T, alpha, Tr, R, m):
#    print ("Z=",Z[0])
    sqrt2 = np.sqrt(2)
    sqrt8 = np.sqrt(8)
    Sdep = R*(-np.log(Z[0]-B)-A/(B*sqrt8)*((m*np.sqrt(Tr))/np.sqrt(alpha))*np.log((Z[0]+(1-sqrt2)*B)/(Z[0]+(1-sqrt2)*B)))
    return Sdep


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"All function to forward calculate and back calculate thermodynamics properties are in this section"

def find_T_rhoP(rho, p): #done
    ## Initial guess for P from ideal gas equation
    R = 8.3145
    mol_mass = 2.01588*10**-3
    T = (p*mol_mass)/(rho*R)
    T1 = T - (tol/100.)*T
    T2 = T + (tol/100.)*T
    T_final = regulaFalsi_T(T1, T2, p, rho)
    return T_final

def find_S_PT(p,T): #done
    H, Href, S, Sref = Cal_Therm(p, T)
    return S, Sref

def find_H_PT(p,T): #done
    H, Href, S, Sref = Cal_Therm(p, T)
    return H, Href

def find_rho_PT(p,T): #done
    mol_mass = 2.01588*10**-3
    R = 8.3145
    Z, phi = PR_H2_Main(T, p)
    rho_val = (p/(R*T))*(mol_mass/Z[0]) 
    return rho_val
    
def find_P_rhoT(rho, T): #done
    ## Initial guess for P from ideal gas equation
    R = 8.3145
    #tol = 15.0
    mol_mass = 2.01588*10**-3
    P = rho*(R/mol_mass)*T ## Initial guess
    P1 = P - (tol/100.)*P
    P2 = P + (tol/100.)*P
    P_final = regulaFalsi(P1, P2, T=T, rho=rho)
    return P_final

def find_T_PS(p, S): #done
    ## Initial guess for P from ideal gas equation
    R = 8.3145
    cp = 11909.0
    #tol = 15.0
    mol_mass = 2.01588*10**-3
    T = Tref*mt.exp((S+R*np.log(p/Pref))/cp) #Initial Guess
    T1 = T - (tol/100.)*T
    T2 = T + (tol/100.)*T
    T_final = regulaFalsi_S(T1, T2, p, S)
    return T_final

def find_S_rhoT(rho, T): #done
    P_cal = find_P_rhoT(rho, T)
    print (P_cal)
    H, Href, S, Sref = Cal_Therm(P_cal, T)
    return S[0]

def find_H_rhoT(rho, T): #done
    P_cal = find_P_rhoT(rho, T)
    H, Href, S, Sref = Cal_Therm(P_cal, T)
    return H[0]

def find_rho_PS(p, s): #done
    T_cal = find_T_PS(p, s)
    mol_mass = 2.01588*10**-3
    R = 8.3145
    Z, phi = PR_H2_Main(T_cal, p)
    rho_cal = (p/(R*T_cal))*(mol_mass/Z[0]) 
    return rho_cal
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

" All function to calculate objective function for Back Calculations of Thermodynamics prop. are here in this section"    

def func(P, T, rho):
    ## for P
    R = 8.3145
    mol_mass = 2.01588*10**-3
    Z, phi = PR_H2_Main(T, P)
    X = P/Z - rho*R*T/mol_mass
    ## for T
    return X


def func_TPrho(T, p, rho):
    ## for P
    R = 8.3145
    mol_mass = 2.01588*10**-3
    Z, phi = PR_H2_Main(T, p)
    X = p/Z[0] - rho*R*T/mol_mass
    ## for T
    return X


def func_H_S(T, P, S):
    ## for P
    R = 8.3145
    mol_mass = 2.01588*10**-3
    A, B, alpha, Tr, m = PR_H2_Coeff(T, P)
    Z, phi = PR_H2_Main(T, P)
    H_g, H_ref, S_g, S_ref= Cal_Therm(P, T)
    X = S - S_g
    return X

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

" All funnction to to calculate roots of objective funstion with false root methods are here in thsi section"
    
def regulaFalsi( a , b, T, rho ):
    MAX_ITER = 1000000
    if func(a, T, rho) * func(b, T, rho) >= 0:
        print("You have not assumed right a and b")
        return -1
     
    c = a # Initialize result
     
    for i in range(MAX_ITER):
        # Find the point that touches x axis
        c = (a * func(b, T, rho) - b * func(a, T, rho))/ (func(b, T, rho) - func(a, T, rho))
         
        # Check if the above found point is root
        if func(c, T, rho) == 0:
            break
         
        # Decide the side to repeat the steps
        elif func(c, T, rho) * func(a, T, rho) < 0:
            b = c
        else:
            a = c
        return c
    
    
def regulaFalsi_T( a , b, p, rho ):
    MAX_ITER = 1000000
    if func_TPrho(a, p, rho) * func_TPrho(b, p, rho) >= 0:
        print("You have not assumed right a and b")
        return -1
     
    c = a # Initialize result
     
    for i in range(MAX_ITER):
        # Find the point that touches x axis
        c = (a * func_TPrho(b, p, rho) - b * func_TPrho(a, p, rho))/ (func_TPrho(b, p, rho) - func_TPrho(a, p, rho))
         
        # Check if the above found point is root
        if func_TPrho(c, p, rho) == 0:
            break
         
        # Decide the side to repeat the steps
        elif func_TPrho(c, p, rho) * func_TPrho(a, p, rho) < 0:
            b = c
        else:
            a = c
        return c

    
    
def regulaFalsi_S( a , b, p, s):
    MAX_ITER = 1000000
    if func_H_S(a, p, s) * func_H_S(b, p, s) >= 0:
        print("You have not assumed right a and b")
        return -1
     
    c = a # Initialize result
     
    for i in range(MAX_ITER):
        # Find the point that touches x axis
        c = (a * func_H_S(b, p, s) - b * func_H_S(a, p, s))/ (func_H_S(b, p, s) - func_H_S(a, p, s))
         
        # Check if the above found point is root
        if func_H_S(c, p, s) == 0:
            break
         
        # Decide the side to repeat the steps
        elif func_H_S(c, p, s) * func_H_S(a, p, s) < 0:
            b = c
        else:
            a = c
        return c
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++