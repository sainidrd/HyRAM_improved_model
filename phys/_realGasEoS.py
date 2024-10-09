# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 12:27:40 2022

@author: saini
"""

class PengRobinsonH2N2:
    def __init__(self, tol=35.0):
        self.tol = tol
        self.R = 8.3145
        self.Tcrit = np.array([33.2, 126.2, 154.581])  # Critical temperature in K
        self.Pcrit = np.array([1.297e6, 3.394e6, 5.043e6])  # Critical pressure in Pa
        self.omega = np.array([-0.219, 0.039, 0.0222])
        
    def PR_H2_N2_Coeff(T, P, x_H2):
        x_Air = 1-x_H2
        x_N2 = x_Air 
        x_O2 = 0.0 
        X = np.array([x_H2, x_N2, x_O2])
        Tcrit = np.array([33.2, 126.2, 154.581]) #in kelvin
        Pcrit = np.array([1.297*10**6, 3.394*10**6, 5.0430*10**6]) # in Pa
        omega = np.array([-0.219,  0.039, 0.0222])
        R = 8.3145
        a = np.zeros(3)
        b = np.zeros(3)
        Tr = np.zeros(3)
        m = np.zeros(3)
        alpha = np.zeros(3)
        for i in range (3):
            Tr[i]    =  T/Tcrit[i]
            m[i]     =  0.37464 + 1.54226*omega[i] - 0.26992*omega[i]**2
            alpha[i] =  (1+m[i]*(1-np.sqrt(Tr[i])))**2
            a[i]     =  0.45724*alpha[i]*((R*Tcrit[i])**2)/Pcrit[i]
            b[i]     =  0.07780*(R*Tcrit[i])/Pcrit[i]
        Ka_New =   [[0.0, 0.0436, 0.433],
                    [0.0436, 0.0, -0.0821],
                    [0.433, -0.0821, 0.0]]
        Kb_New =   [[0.0, -0.1611, 0.0382],
                    [-0.1611, 0.0, -0.2633],
                    [0.0382, -0.2633, 0.0]]
        a_m = np.zeros((3,3))
        b_m = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                a_m[i][j]  =  (1-Ka_New[i][j])*(X[i]*X[j])*(a[i]*a[j])**0.5
                b_m[i][j]  =  (1-Kb_New[i][j])*(X[i]*X[j])*((b[i]+b[j])/2)
        A     =  (np.sum(a_m)*P)/(R*T)**2
        B     =  (np.sum(b_m)*P)/(R*T)
        return A, B, alpha, Tr, m
    
    def PR_H2_N2_Coeff_ab(T, P, x_H2):
        x_Air = 1-x_H2
        x_N2 = x_Air #(0.767)*x_Air
        x_O2 = 0.0 #(0.233)*x_Air
        X = np.array([x_H2, x_N2, x_O2])
        Tcrit = np.array([33.2, 126.2, 154.581]) #in kelvin
        Pcrit = np.array([1.297*10**6, 3.394*10**6, 5.0430*10**6]) # in Pa
        omega = np.array([-0.219,  0.039, 0.0222])
        R = 8.3145
        a = np.zeros(3)
        b = np.zeros(3)
        Tr = np.zeros(3)
        m = np.zeros(3)
        alpha = np.zeros(3)
        for i in range (3):
            Tr[i] = T/Tcrit[i]
            m[i] = 0.37464 + 1.54226*omega[i] - 0.26992*omega[i]**2
            alpha[i] = (1+m[i]*(1-np.sqrt(Tr[i])))**2
            a[i] =  0.45724*alpha[i]*((R*Tcrit[i])**2)/Pcrit[i]
            b[i]     = 0.07780*(R*Tcrit[i])/Pcrit[i]
        Ka_New =   [[0.0, 0.0436, 0.433],
                    [0.0436, 0.0, -0.0821],
                    [0.433, -0.0821, 0.0]]
        Kb_New =   [[0.0, -0.1611, 0.0382],
                    [-0.1611, 0.0, -0.2633],
                    [0.0382, -0.2633, 0.0]]
        a_m = np.zeros((3,3))
        b_m = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                a_m[i][j]  =  (1-Ka_New[i][j])*(X[i]*X[j])*(a[i]*a[j])**0.5
                b_m[i][j]  =  (1-Kb_New[i][j])*(X[i]*X[j])*((b[i]+b[j])/2)
        
        return np.sum(a_m), np.sum(b_m)
    
    def PR_H2_N2_Coeff_ab_ind(T, P, x_H2):
        x_Air = 1-x_H2
        x_N2 = x_Air #(0.767)*x_Air
        x_O2 = 0.0 #(0.233)*x_Air
        X = np.array([x_H2, x_N2, x_O2])
        Tcrit = np.array([33.2, 126.2, 154.581]) #in kelvin
        Pcrit = np.array([1.297*10**6, 3.394*10**6, 5.0430*10**6]) # in Pa
        omega = np.array([-0.219,  0.039, 0.0222])
        R = 8.3145
        a = np.zeros(3)
        b = np.zeros(3)
        Tr = np.zeros(3)
        m = np.zeros(3)
        alpha = np.zeros(3)
        for i in range (3):
            Tr[i]    =  T/Tcrit[i]
            m[i]     =  0.37464 + 1.54226*omega[i] - 0.26992*omega[i]**2
            alpha[i] =  (1+m[i]*(1-np.sqrt(Tr[i])))**2
            a[i]     =  0.45724*alpha[i]*((R*Tcrit[i])**2)/Pcrit[i]
            b[i]     =  0.07780*(R*Tcrit[i])/Pcrit[i]
        return a, b
    
    def PR_H2_N2_Coeff_deriv(T, P, x_H2):
        # Constants and inputs
        Tcrit = np.array([33.2, 126.2])  # in Kelvin
        Pcrit = np.array([1.297e6, 3.394e6])  # in Pa
        omega = np.array([-0.219, 0.039])
        R = 8.3145
        
        # Mixture composition
        x_N2 = 1 - x_H2
        
        # Redlich-Kwong parameters
        Tr = T / Tcrit
        alpha = (1 + 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2 * (1 - np.sqrt(Tr))) ** 2
        a = 0.45724 * alpha * (R ** 2 * Tcrit ** 2) / Pcrit
        b = 0.07780 * R * Tcrit / Pcrit
        
        # Calculate compressibility factor Z
        A = a[0] * x_H2 ** 2 + a[1] * x_N2 ** 2 + 2 * x_H2 * x_N2 * np.sqrt(a[0] * a[1])
        B = b[0] * x_H2 + b[1] * x_N2
        #Z = P * (R * T) / (Z * (R * T) * (Z - B)) - A / (Z * (Z + B) * (Z + B))
        Z, phi = PR_H2_N2_Main(T, P, x_H2)
        # Calculate derivatives
        dA_dT = 0.45724 * (1 + 1.54226 * omega - 0.26992 * omega ** 2) * (R ** 2 * Tcrit ** 2) / (Pcrit * Tcrit) * \
                (-0.37464 * (1 - np.sqrt(Tr)) - 0.5 * 1.54226 * omega / np.sqrt(Tr) + 2 * 0.26992 * omega * np.sqrt(Tr))
        dB_dT = 0.07780 * (R * Tcrit) / Pcrit
        dA_dP = -0.45724 * (1 + 1.54226 * omega - 0.26992 * omega ** 2) * (R * Tcrit) ** 2 / (Pcrit ** 2) / Tcrit
        dB_dP = 0.07780 * (R * Tcrit) / (Pcrit ** 2)
        
        return A, B, dA_dT, dB_dT, dA_dP, dB_dP
    
    
    def find_real_roots(r):
        N = len(r)
        rr=[]
        RR = np.zeros((N))
        for i in range(N):
            if r[i].imag==0.0:        
                if r[i].real > 0.0:
                   RR[i] = r[i]
        for i in range(N):
            if RR[i] > 0.0:
               r_positive = RR[i]
               rr.append(r_positive)
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
            #print ("More than one roots exist for PR equations for H2-N2 Mixtures")
            RR=[]
            RR.append(real_roots[0])
            #print ("RR=", RR)
            phi =  fugacity(RR, A, B)
            real_roots = RR
    #        XX  = np.where(real_roots > B )[0] # index where roots are > B
    #        real_RR=[]
    #        for i in range(len(XX)):
    #         #   print ("real_roots[XX]", real_roots[XX[i]])
    #            ROOt = real_roots[XX[i]]
    #            real_RR.append(ROOt)
    #        #print("ROOt", real_RR)
    #        phi = fugacity(real_RR, A, B) # find fugacity of those roots
    #        print ("phi=",phi)
    #        if len(phi) > 1:
    #            XXX = np.where(phi == min(phi))[0]
    #         #   print ("min fug", XXX[0])
    #            real_roots = real_RR[XXX[0]]
    #          #  print ("real_root min fug",real_roots)
    #            phi = phi[XXX[0]]
        elif len(real_roots)==1:
            #print ("Just one root exist for PR equation for H2-N2 Mixtures")
            phi = fugacity(real_roots, A, B)
        return real_roots, phi
    
    def PR_H2_N2_Main(T, P, x_H2):
        
        #Calculate Coeffs for Peng-Robinson Polynomial
        A, B, alpha, Tr, m = PR_H2_N2_Coeff(T, P, x_H2 )
        
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
        #print ("real_roots are",real_roots)
        ZZ , phi_ZZ = Sort_roots(real_roots, A, B)
        #print("After roots sorting ZZ=", ZZ)
        return ZZ, phi_ZZ
    
    def Cal_Therm(P, T, xH2):
        Tcrit = np.array([33.2, 126.2]) #in kelvin
        Pcrit = np.array([1.297*10**6, 3.394*10**6]) # in Pa
        Tref = xH2*Tcrit[0]+(1-xH2)*Tcrit[1]
        Pref = xH2*Pcrit[0]+(1-xH2)*Pcrit[1]
        xN2 = 1-xH2
        R = 8.3145 #2.016
        Cp_Coeff_N2 = np.array([-2.07E-09/6, 2.19E-06/5, -9.09E-04/4, 1.86E-01/3, -1.89E+01/2, 1.81E+03])
        Cp_Coeff_H2 = np.array([-2.64032606e-15/9,  4.34339102e-12/8, -2.94216657e-09/7, \
                                1.04846126e-06/6, -2.04658359e-04/5,  1.94881968e-02/4, \
                                -3.77290489e-01/3, -3.89881763e+01/2, 1.19295272e+04])
        Cp_Coeff_O2 = np.array([-7.74E-08/6, 7.91E-05/5, -3.14E-02/4, 6.06E+00/3, -5.66E+02/2, 21311.8229])
        H_ideal = (Cp_Coeff_N2[0]*(T**6- Tref**6) + Cp_Coeff_N2[1]*(T**5- Tref**5) + Cp_Coeff_N2[2]*(T**4- Tref**4) + \
                  Cp_Coeff_N2[3]*(T**3- Tref**3) + Cp_Coeff_N2[4]*(T**2- Tref**2) + Cp_Coeff_N2[5]*(T- Tref))*xN2  + \
                  (Cp_Coeff_H2[0]*(T**9- Tref**9) + Cp_Coeff_H2[1]*(T**8- Tref**8) + Cp_Coeff_H2[2]*(T**7- Tref**7) + \
                   Cp_Coeff_H2[3]*(T**6- Tref**6) + Cp_Coeff_H2[4]*(T**5- Tref**5) + Cp_Coeff_H2[5]*(T**4- Tref**4) + \
                   Cp_Coeff_H2[6]*(T**3- Tref**3) + Cp_Coeff_H2[7]*(T**2- Tref**2) + Cp_Coeff_H2[8]*(T- Tref))*xH2 
        # Departure from ideal gas behaviour
        Z, phi = PR_H2_N2_Main(T, P, xH2)
        print ("Z for mixtures=", Z)
        A, B, alpha, Tr, m = PR_H2_N2_Coeff(T, P, xH2)
        H_dep = Cal_H_dep(Z, A, B, T, alpha, Tr, R, m)
        print ("Hdep=", H_dep)
        print ("Hideal=", H_ideal)
        Z_ref, phi_ref = PR_H2_N2_Main(Tref, Pref, xH2)
        H_ref = Cal_H_dep(Z_ref, A, B, Tref, alpha, Tr, R, m)
        print ("Href=", H_ref)
        H_dep = xH2*H_dep[0] + xN2*H_dep[1]
        H_ref = xH2*H_ref[0] + xN2*H_ref[1]
        print ("H_dep final=", H_dep)
        print ("H_ref final=", H_ref)
        H = H_dep + H_ideal - H_ref
        print ("H_final=", H)
        return H, H_ref
    
    
    def Cal_Therm_new(P, T, xH2):
        Tcrit = np.array([33.2, 126.2]) #in kelvin
        Pcrit = np.array([1.297*10**6, 3.394*10**6]) # in Pa
        Tref = xH2*Tcrit[0]+(1-xH2)*Tcrit[1]
        Pref = xH2*Pcrit[0]+(1-xH2)*Pcrit[1]
        xN2 = 1-xH2
        R = 8.3145 #2.016
        Cp_Coeff_N2 = np.array([-2.07E-09/6, 2.19E-06/5, -9.09E-04/4, 1.86E-01/3, -1.89E+01/2, 1.81E+03])
        Cp_Coeff_H2 = np.array([-2.64032606e-15/9,  4.34339102e-12/8, -2.94216657e-09/7, \
                                1.04846126e-06/6, -2.04658359e-04/5,  1.94881968e-02/4, \
                                -3.77290489e-01/3, -3.89881763e+01/2, 1.19295272e+04])
        Cp_Coeff_O2 = np.array([-7.74E-08/6, 7.91E-05/5, -3.14E-02/4, 6.06E+00/3, -5.66E+02/2, 21311.8229])
        H_ideal = (Cp_Coeff_N2[0]*(T**6- Tref**6) + Cp_Coeff_N2[1]*(T**5- Tref**5) + Cp_Coeff_N2[2]*(T**4- Tref**4) + \
                  Cp_Coeff_N2[3]*(T**3- Tref**3) + Cp_Coeff_N2[4]*(T**2- Tref**2) + Cp_Coeff_N2[5]*(T- Tref))*xN2  + \
                  (Cp_Coeff_H2[0]*(T**9- Tref**9) + Cp_Coeff_H2[1]*(T**8- Tref**8) + Cp_Coeff_H2[2]*(T**7- Tref**7) + \
                   Cp_Coeff_H2[3]*(T**6- Tref**6) + Cp_Coeff_H2[4]*(T**5- Tref**5) + Cp_Coeff_H2[5]*(T**4- Tref**4) + \
                   Cp_Coeff_H2[6]*(T**3- Tref**3) + Cp_Coeff_H2[7]*(T**2- Tref**2) + Cp_Coeff_H2[8]*(T- Tref))*xH2 
        # Departure from ideal gas behaviour
        Z, phi = PR_H2_N2_Main(T, P, xH2)
        print ("Z for mixtures=", Z)
        A, B, alpha, Tr, m = PR_H2_N2_Coeff(T, P, xH2)
        H_dep = Cal_H_dep_new(Z, A, B, T, alpha, Tr, R, m, xH2, P)
        print ("Hdep=", H_dep)
        print ("Hideal=", H_ideal)
        Z_ref, phi_ref = PR_H2_N2_Main(Tref, Pref, xH2)
        H_ref = Cal_H_dep_new(Z_ref, A, B, Tref, alpha, Tr, R, m , xH2, P)
        print ("Href=", H_ref)
        #H_dep = xH2*H_dep[0] + xN2*H_dep[1]
        #H_ref = xH2*H_ref[0] + xN2*H_ref[1]
        print ("H_dep final=", H_dep)
        print ("H_ref final=", H_ref)
        H = H_dep + H_ideal - H_ref
        print ("H_final=", H)
        return H_dep, H_ref
    
    def Cal_H_dep(Z, A, B, T, alpha, Tr, R, m):
        sqrt2 = np.sqrt(2)
        print ("m=", m)
        print ("alpha=", alpha)
        print ("Tr=", Tr)
        try:
            Hdep = R*T*(Z-1.-A/B/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))* \
                        (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
        except:
            Hdep = R*T*(Z[0]-1.-A/B/2/sqrt2*np.log((Z[0]+(1+sqrt2)*B)/(Z[0]+(1-sqrt2)*B))* \
                        (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
            #Hdep = R*T*(Z-1.-A/B/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))* \
            #            (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
        return Hdep
    
    
    def Cal_H_dep_new(Z, A, B, T, alpha, Tr, R, m, xH2, P):
        x_H2 = xH2
        x_N2 = 1-x_H2
        Tcrit = np.array([33.2, 126.2]) #in kelvin
        Pcrit = np.array([1.297*10**6, 3.394*10**6]) # in Pa
        Tr    = T/Tcrit
        omega = np.array([-0.219,  0.039])
        R = 8.3145
        m =  0.37464 + 1.54226*omega - 0.26992*omega**2
        alpha = (1+m*(1-np.sqrt(Tr)))**2
        a     = 0.45724*alpha*((R*Tcrit)**2)/Pcrit
        b     = 0.07780*(R*Tcrit)/Pcrit
        k_12  = -0.0878 #from ref 101
        k_11  = 0.0
        k_22  = 0.0
        a11   = a[0]
        a22   = a[1]
        a12   = (1.-k_12)*(a[0]*a[1])**0.5
        b11   = b[0]
        b22   = b[1]
        b12   = (b[0]+b[1])/2.0 # From Clapeyron.jl code
        a_m   = x_H2*x_H2*a11 + x_N2*x_N2*a22 + 2.0*x_H2*x_N2*a12
        b_m   = x_H2*x_H2*b11 + x_N2*x_N2*b22 + x_H2*x_N2*b12
        A     = (a_m*P)/(R*T)**2
        B     = (b_m*P)/(R*T)
        sqrt2 = np.sqrt(2)
        print ("a_m=", a_m)
        print ("alpha=", alpha)
        print ("Z=", Z[0])
        Hdep = R*T*(Z[0]-1.-1.0/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))*\
                    -(x_H2*x_H2*a[0]+x_N2*x_N2*a[1]+2.0*x_H2*x_N2*(1-k_12)*(a[0]*a[1])**0.5)-\
                        ((x_H2*x_H2*(a[0]/(Tcrit[0]*Tr[0]**0.5))*(alpha[0])**0.5)+\
                        (x_N2*x_N2*(a[1]/(Tcrit[1]*Tr[1]**0.5))*(alpha[1])**0.5)+\
                        (x_H2*x_N2*(1-k_12)*\
                        ((a[0]**0.5/a[1]**0.5)*a[1]*omega[1]*(alpha[1]**0.5)/Tcrit[1]*(Tr[1]**0.5)+
                         (a[1]**0.5/a[0]**0.5)*a[0]*omega[0]*(alpha[0]**0.5)/Tcrit[0]*(Tr[0]**0.5)))\
                                ))
        #try:
        #    Hdep = R*T*(Z-1.-A/B/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))* \
        #                (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
        #except:
        #    Hdep = R*T*(Z[0]-1.-A/B/2/sqrt2*np.log((Z[0]+(1+sqrt2)*B)/(Z[0]+(1-sqrt2)*B))* \
        #                (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
            #Hdep = R*T*(Z-1.-A/B/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))* \
            #            (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
        return Hdep
    
    def Cal_Cp_mixture(T, xH2):
        Cp_Coeff_N2=np.array([1.50848431e-15, -2.46381879e-12,  1.73643875e-09, -6.89801002e-07, \
                              1.69042289e-04, -2.62077000e-02,  2.51795719e+00, -1.37898047e+02, \
                              4.37824054e+03])
        Cp_Coeff_H2 = np.array([-2.64032606e-15,  4.34339102e-12, -2.94216657e-09, \
                                1.04846126e-06, -2.04658359e-04,  1.94881968e-02, \
                                -3.77290489e-01, -3.89881763e+01, 1.19295272e+04])
        Cp_H2 = Cp_Coeff_H2[0]*(T**8) + Cp_Coeff_H2[1]*(T**7) + Cp_Coeff_H2[2]*(T**6) + \
                    Cp_Coeff_H2[3]*(T**5) + Cp_Coeff_H2[4]*(T**4) + Cp_Coeff_H2[5]*(T**3) + \
                    Cp_Coeff_H2[6]*(T**2) + Cp_Coeff_H2[7]*(T**1) + Cp_Coeff_H2[8]
        Cp_N2 = (Cp_Coeff_N2[0]*(T**8) + Cp_Coeff_N2[1]*(T**7) + Cp_Coeff_N2[2]*(T**6) + \
                 Cp_Coeff_N2[3]*(T**5) + Cp_Coeff_N2[4]*(T**4) + Cp_Coeff_N2[5]*(T**3) + \
                 Cp_Coeff_N2[6]*(T**2) + Cp_Coeff_N2[7]*(T**1) + Cp_Coeff_N2[8])
        Cp_mix = Cp_H2*xH2 + (1-xH2)*Cp_N2
        return Cp_mix
    
    
    def PR_H2_N2_Derivatives(T, P, x_H2):
        x_N2 = 1 - x_H2
        Tcrit = np.array([33.2, 126.2])  # in Kelvin
        Pcrit = np.array([1.297e6, 3.394e6])  # in Pa
        Tr = T / Tcrit
        omega = np.array([-0.219, 0.039])
        R = 8.3145
        
        # Calculate derivatives of A and B with respect to temperature and pressure
        dA_dT = 0.45724 * (1 + 1.54226 * omega - 0.26992 * omega ** 2) * (R ** 2 * Tcrit ** 2) / (Pcrit * Tcrit) * \
                (-0.37464 * (1 - np.sqrt(Tr)) - 0.5 * 1.54226 * omega / np.sqrt(Tr) + 2 * 0.26992 * omega * np.sqrt(Tr))
        dB_dT = 0.07780 * (R * Tcrit) / Pcrit
        dA_dP = -0.45724 * (1 + 1.54226 * omega - 0.26992 * omega ** 2) * (R * Tcrit) ** 2 / (Pcrit ** 2) / Tcrit
        dB_dP = 0.07780 * (R * Tcrit) / (Pcrit ** 2)
        
        return dA_dT, dB_dT, dA_dP, dB_dP
    
    def find_rho_PT_Mixture(p,T, x_H2): #done
        mol_mass = ((x_H2*2.016 + (1.0-x_H2)*28.010))*10**-3
        R = 8.3145
        Z, phi = PR_H2_N2_Main(T, p, x_H2)
        try:
            rho_val = (p/(R*T))*(mol_mass/Z) 
        except:
            rho_val = (p/(R*T))*(mol_mass/Z[0]) 
        return rho_val
    
    def find_H_PT(p,T, x_H2): #done
        H, Href = Cal_Therm(p, T, x_H2)
        return H
    
    def find_H_PT_new(p,T, x_H2): #done
        H, Href = Cal_Therm_new(p, T, x_H2)
        return H[0]
    
    
    def find_T_rhoP_Mixture(rho, p, xH2): #done
        ## Initial guess for P from ideal gas equation
        R = 8.314462618 #8.3145
        tol = 10.0
        MW_H2 = 2.01588*10**-3
        MW_N2 = 0.02896546 #28.010*10**-3
        Y_ = xH2*MW_H2/(xH2*MW_H2+(1-xH2)*MW_N2)
        mol_mass = MW_N2*MW_H2/(Y_*(MW_N2-MW_H2) + MW_H2)
        #mol_mass = 1/(xH2/MW_H2 + (1-xH2)/MW_N2)
        T = (p*mol_mass)/(rho*R)
        #print ("T=", T)
        #T_cl0 = ambient.P*MW_cl0/(const.R*self.initial_node.rho_cl)
        T1 = T - (tol/100.)*T
        T2 = T + (tol/100.)*T
        #print ("T1=", T1)
        #print ("T2=", T2)
        #print ("xH2=", xH2)
        T_final = regulaFalsi_T(T1, T2, p, rho, xH2)
        return T_final
    
    def derivative_T_rhoP_Mixture(rho, p, xH2):
        # Calculate T using existing function
        MW_H2 = 2.01588*10**-3
        MW_N2 = 0.02896546
        P = 101325.0
        R = 8.314462618 #8.3145
        Y_ = xH2*MW_H2/(xH2*MW_H2+(1-xH2)*MW_N2)
        mol = MW_N2*MW_H2/(Y_*(MW_N2-MW_H2) + MW_H2)
        T = find_T_rhoP_Mixture(rho, p, xH2)
        a, b = PR_H2_N2_Coeff_ab(T, P, xH2)
        M = mol
        dT_drho = 2*a*(mol/R)*((mol**2 - (b*rho)**2)/(mol**2 + 2*mol*b*rho - (b*rho)**2)**2)\
                  - ((mol*T)/rho)*(1./(mol-b*rho)) 
        return dT_drho
    
    def derivative_T_Y_Mixture(rho, p, xH2):
        # Calculate T using existing function
        MW_H2 = 2.01588*10**-3
        MW_N2 = 0.02896546 
        R = 8.314462618 #8.3145
        P = 101325.0
        Y_ = xH2*MW_H2/(xH2*MW_H2+(1-xH2)*MW_N2)
        mol = MW_N2*MW_H2/(Y_*(MW_N2-MW_H2) + MW_H2)
        T = find_T_rhoP_Mixture(rho, p, xH2)
        dmoldY = ((MW_H2-MW_N2)*(MW_H2*MW_N2))/(MW_H2 + Y_*(MW_N2-MW_H2))**2
        am, bm = PR_H2_N2_Coeff_ab(T, P, xH2)
        a, b = PR_H2_N2_Coeff_ab_ind(T, P, xH2)
        F_N2 = 0.767
        F_O2 = 0.237
        Ka_New =   [[0.0, 0.0436, 0.433],
                    [0.0436, 0.0, -0.0821],
                    [0.433, -0.0821, 0.0]]
        Kb_New =   [[0.0, -0.1611, 0.0382],
                    [-0.1611, 0.0, -0.2633],
                    [0.0382, -0.2633, 0.0]]
        dadY = (2*xH2*a[0]*(1-Ka_New[0][0]) - 2*F_N2*(1-xH2)*a[1]*(1-Ka_New[1][1]) - 2*a[2]*(F_O2**2)*(1-xH2)*(1-Ka_New[2][2])\
                + 2*(1-Ka_New[0][1])*F_N2*(1-2*xH2)*(a[0]*a[1])**0.5 + 2*(1-Ka_New[0][2])*F_O2*(1-2*xH2)*(a[0]*a[2])**0.5 \
                - 4*(1-Ka_New[1][2])*F_N2*F_O2*(1-xH2)*(a[1]*a[2])**0.5)*((MW_H2*MW_N2)/(MW_H2 + Y_*(MW_N2-MW_H2))**2)
        dbdY = (2*xH2*b[0]*(1-Kb_New[0][0]) - 2*F_N2*(1-xH2)*b[1]*(1-Kb_New[1][1]) - 2*F_O2*(1-xH2)*(1-Kb_New[2][2])*b[2] \
                + (1-Kb_New[0][1])*F_N2*(1-2*xH2)*(b[0]+b[1]) + (1-Kb_New[0][2])*F_O2*(1-2*xH2)*(b[0]+b[2])\
                    - (1-Kb_New[1][2])*F_N2*F_O2*(1-xH2)*(b[1]+b[2]))*((MW_H2*MW_N2)/(MW_H2 + Y_*(MW_N2-MW_H2))**2)
                
        #        - 2*(1-xH2)*b22 + b12*(1-2*xH2))*((MW_H2*MW_N2)/(MW_H2 + Y_*(MW_N2-MW_H2))**2)
        dTdY = (T/(mol-bm*rho))*(dmoldY - rho*dbdY)\
            + ((rho/R)*((mol-bm*rho)/(mol**2 + 2*mol*bm*rho - (bm*rho)**2)))*(dadY - ((2*am)/(mol**2 + 2*mol*bm*rho - (bm*rho)**2))\
                                                                           *((mol+bm*rho)*dmoldY + (mol-bm*rho)*rho*dbdY))
          #  ((a*rho**2)/((mol**2 + 2*mol*b*rho - (b*rho)**2)**2))*((mol+b*rho)*2*dmoldY + (rho*mol-b*rho**2)*2*dbdY)\
          #  - (dadY*rho**2)/((mol**2 + 2*mol*b*rho - (b*rho)**2))
        return dTdY
    
    
    def derivative_T_rhoP_Mixture_H(rho, p, xH2):
        # Calculate T using existing function
        T = find_T_rhoP_Mixture(rho, p, xH2)
        
        # Perturb rho slightly for finite difference calculation
        delta_rho = 1e-8 * rho  # Smaller perturbation for higher accuracy
        rho_plus = rho + delta_rho
        rho_minus = rho - delta_rho
        
        # Calculate T for rho+ and rho- using existing function
        T_plus = find_T_rhoP_Mixture(rho_plus, p, xH2)
        T_minus = find_T_rhoP_Mixture(rho_minus, p, xH2)
        
        # Calculate derivative using central difference approximation
        dT_drho = (-T_plus + 8 * T - 8 * T_minus + T_plus) / (12 * delta_rho)
        
        return dT_drho
    
    def dT_drho_PR_H2_N2(T, P, x_H2):
        R = 8.314462618 #8.3145
        # Calculate derivatives of A and B with respect to temperature and pressure
        #dA_dT, dB_dT, dA_dP, dB_dP = PR_H2_N2_Derivatives(T, P, x_H2)
        A, B, dA_dT, dB_dT, dA_dP, dB_dP = PR_H2_N2_Coeff_deriv(T, P, x_H2)
        #dA_dT = dA_dT
        #dB_dT = dB_dT
        #dA_dP = dA_dP
        #dB_dP = dB_dP
        
        # Calculate derivatives of Z with respect to temperature and pressure
        dZ_dT = -(dA_dT - 2*P*dA_dP)/(R*T**2)
        dZ_dP = (B - 2*P*dB_dP)/(R*T)
        Z, phi = PR_H2_N2_Main(T, P, x_H2)
        # Calculate derivatives of T with respect to density using chain rule
        #dT_drho = dT_dZ*dZ_dT + dT_dP*dP_dZ*dZ_dT + dT_dP*dP_dZ*dZ_dP
        try:
            #rho_val = (p/(R*T))*(mol_mass/Z) 
            dT_drho = -T * dZ_dT / (Z * (Z + B))
        except:
            #rho_val = (p/(R*T))*(mol_mass/Z[0])
            dT_drho = -T * dZ_dT / (Z[0] * (Z[0] + B))
        ##dT_drho = -T * dZ_dT / (Z * (Z + B))
        return dT_drho
    
    def find_Cp_H2_N2_Mixture(T, xH2):
        Cp_mix = Cal_Cp_mixture(T, xH2)
        return Cp_mix
    
    def find_T_PH_Mixture(p, H, xH2): #done
        ## Initial guess for P from ideal gas equation
        R = 8.3145
        cp_N2 = 4.37824054e+03
        cp_H2 = 1.19295272e+04
        Tcrit = np.array([33.2, 126.2]) #in kelvin
        Pcrit = np.array([1.297*10**6, 3.394*10**6]) # in Pa
        Tref = xH2*Tcrit[0]+(1-xH2)*Tcrit[1]
        Pref = xH2*Pcrit[0]+(1-xH2)*Pcrit[1]
        #Tref =  33.2*xH2 + 126.2*(1-xH2)
        cp = xH2*cp_H2 + (1-xH2)*cp_N2
        T = H/cp + Tref
        #print ("Assumed values of T=", T)
        tol = 25
        T1 = T - (tol/100.)*T
        T2 = T + (tol/100.)*T
        #print ("T_1=", T1)
        #print ("T_2=", T2)
        T_final = regulaFalsi_H(T1, T2, p, H, xH2)
        print ("T_final=", T_final)
        return T_final
    
    
    def func_TPrho(T, p, rho, xH2):
        ## for P
        R = 8.3145
        MW_H2 = 2.01588*10**-3
        MW_N2 = 0.02896546 #28.010*10**-3
        Y_ = xH2*MW_H2/(xH2*MW_H2+(1-xH2)*MW_N2)
        mol_mass = MW_N2*MW_H2/(Y_*(MW_N2-MW_H2) + MW_H2)
        #mol_mass = 1/(xH2/MW_H2 + (1-xH2)/MW_N2)
        Z, phi = PR_H2_N2_Main(T, p, xH2)
        #print ("ZZZ", Z)
        try:
            X = p/Z - rho*R*T/mol_mass
        except:
            X = p/Z[0] - rho*R*T/mol_mass
        return X


    def func_H_S(T, P, H, xH2):
        ## for P
        R = 8.3145
        mol_mass = 2.01588*10**-3
        A, B, alpha, Tr, m = PR_H2_N2_Coeff(T, P, xH2)
        Z, phi = PR_H2_N2_Main(T, P, xH2)
        H_g, H_ref= Cal_Therm(P, T, xH2)
        X = H - H_g
        return X
    
    def regulaFalsi_H( a , b, p, h, xH2):
        MAX_ITER = 1000000
        #print ("func_H_S(a, p, h, xH2)=", func_H_S(a, p, h, xH2))
        #print ("func_H_S(b, p, h, xH2)=", func_H_S(b, p, h, xH2))
        if func_H_S(a, p, h, xH2) * func_H_S(b, p, h, xH2) >= 0:
            print("You have not assumed right a and b  in pure H2 regulaFalsi_S")
            return -1
         
        c = a # Initialize result
         
        for i in range(MAX_ITER):
            # Find the point that touches x axis
            c = (a * func_H_S(b, p, h, xH2) - b * func_H_S(a, p, h, xH2))/ (func_H_S(b, p, h, xH2) - func_H_S(a, p, h, xH2))
             
            # Check if the above found point is root
            if func_H_S(c, p, h, xH2) == 0:
                break
             
            # Decide the side to repeat the steps
            elif func_H_S(c, p, h, xH2) * func_H_S(a, p, h, xH2) < 0:
                b = c
            else:
                a = c
            return c
        
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
    
        
        






        
        
    



















import numpy as np
global Tcrit, Pcrit , tol
tol = 35.0


def PR_H2_N2_Coeff(T, P, x_H2):
    x_Air = 1-x_H2
    x_N2 = x_Air #(0.767)*x_Air
    x_O2 = 0.0 #(0.233)*x_Air
    X = np.array([x_H2, x_N2, x_O2])
    Tcrit = np.array([33.2, 126.2, 154.581]) #in kelvin
    Pcrit = np.array([1.297*10**6, 3.394*10**6, 5.0430*10**6]) # in Pa
    omega = np.array([-0.219,  0.039, 0.0222])
    R = 8.3145
    a = np.zeros(3)
    b = np.zeros(3)
    Tr = np.zeros(3)
    m = np.zeros(3)
    alpha = np.zeros(3)
    for i in range (3):
        Tr[i]    =  T/Tcrit[i]
        m[i]     =  0.37464 + 1.54226*omega[i] - 0.26992*omega[i]**2
        alpha[i] =  (1+m[i]*(1-np.sqrt(Tr[i])))**2
        a[i]     =  0.45724*alpha[i]*((R*Tcrit[i])**2)/Pcrit[i]
        b[i]     =  0.07780*(R*Tcrit[i])/Pcrit[i]
    Ka_New =   [[0.0, 0.0436, 0.433],
                [0.0436, 0.0, -0.0821],
                [0.433, -0.0821, 0.0]]
    Kb_New =   [[0.0, -0.1611, 0.0382],
                [-0.1611, 0.0, -0.2633],
                [0.0382, -0.2633, 0.0]]
    a_m = np.zeros((3,3))
    b_m = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            a_m[i][j]  =  (1-Ka_New[i][j])*(X[i]*X[j])*(a[i]*a[j])**0.5
            b_m[i][j]  =  (1-Kb_New[i][j])*(X[i]*X[j])*((b[i]+b[j])/2)
    A     =  (np.sum(a_m)*P)/(R*T)**2
    B     =  (np.sum(b_m)*P)/(R*T)
    return A, B, alpha, Tr, m



def PR_H2_N2_Coeff_ab(T, P, x_H2):
    x_Air = 1-x_H2
    x_N2 = x_Air #(0.767)*x_Air
    x_O2 = 0.0 #(0.233)*x_Air
    X = np.array([x_H2, x_N2, x_O2])
    Tcrit = np.array([33.2, 126.2, 154.581]) #in kelvin
    Pcrit = np.array([1.297*10**6, 3.394*10**6, 5.0430*10**6]) # in Pa
    omega = np.array([-0.219,  0.039, 0.0222])
    R = 8.3145
    a = np.zeros(3)
    b = np.zeros(3)
    Tr = np.zeros(3)
    m = np.zeros(3)
    alpha = np.zeros(3)
    for i in range (3):
        Tr[i] = T/Tcrit[i]
        m[i] = 0.37464 + 1.54226*omega[i] - 0.26992*omega[i]**2
        alpha[i] = (1+m[i]*(1-np.sqrt(Tr[i])))**2
        a[i] =  0.45724*alpha[i]*((R*Tcrit[i])**2)/Pcrit[i]
        b[i]     = 0.07780*(R*Tcrit[i])/Pcrit[i]
    Ka_New =   [[0.0, 0.0436, 0.433],
                [0.0436, 0.0, -0.0821],
                [0.433, -0.0821, 0.0]]
    Kb_New =   [[0.0, -0.1611, 0.0382],
                [-0.1611, 0.0, -0.2633],
                [0.0382, -0.2633, 0.0]]
    a_m = np.zeros((3,3))
    b_m = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            a_m[i][j]  =  (1-Ka_New[i][j])*(X[i]*X[j])*(a[i]*a[j])**0.5
            b_m[i][j]  =  (1-Kb_New[i][j])*(X[i]*X[j])*((b[i]+b[j])/2)
    
    return np.sum(a_m), np.sum(b_m)



def PR_H2_N2_Coeff_ab_ind(T, P, x_H2):
    x_Air = 1-x_H2
    x_N2 = x_Air #(0.767)*x_Air
    x_O2 = 0.0 #(0.233)*x_Air
    X = np.array([x_H2, x_N2, x_O2])
    Tcrit = np.array([33.2, 126.2, 154.581]) #in kelvin
    Pcrit = np.array([1.297*10**6, 3.394*10**6, 5.0430*10**6]) # in Pa
    omega = np.array([-0.219,  0.039, 0.0222])
    R = 8.3145
    a = np.zeros(3)
    b = np.zeros(3)
    Tr = np.zeros(3)
    m = np.zeros(3)
    alpha = np.zeros(3)
    for i in range (3):
        Tr[i]    =  T/Tcrit[i]
        m[i]     =  0.37464 + 1.54226*omega[i] - 0.26992*omega[i]**2
        alpha[i] =  (1+m[i]*(1-np.sqrt(Tr[i])))**2
        a[i]     =  0.45724*alpha[i]*((R*Tcrit[i])**2)/Pcrit[i]
        b[i]     =  0.07780*(R*Tcrit[i])/Pcrit[i]
    return a, b

def PR_H2_N2_Coeff_deriv(T, P, x_H2):
    # Constants and inputs
    Tcrit = np.array([33.2, 126.2])  # in Kelvin
    Pcrit = np.array([1.297e6, 3.394e6])  # in Pa
    omega = np.array([-0.219, 0.039])
    R = 8.3145
    
    # Mixture composition
    x_N2 = 1 - x_H2
    
    # Redlich-Kwong parameters
    Tr = T / Tcrit
    alpha = (1 + 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2 * (1 - np.sqrt(Tr))) ** 2
    a = 0.45724 * alpha * (R ** 2 * Tcrit ** 2) / Pcrit
    b = 0.07780 * R * Tcrit / Pcrit
    
    # Calculate compressibility factor Z
    A = a[0] * x_H2 ** 2 + a[1] * x_N2 ** 2 + 2 * x_H2 * x_N2 * np.sqrt(a[0] * a[1])
    B = b[0] * x_H2 + b[1] * x_N2
    #Z = P * (R * T) / (Z * (R * T) * (Z - B)) - A / (Z * (Z + B) * (Z + B))
    Z, phi = PR_H2_N2_Main(T, P, x_H2)
    # Calculate derivatives
    dA_dT = 0.45724 * (1 + 1.54226 * omega - 0.26992 * omega ** 2) * (R ** 2 * Tcrit ** 2) / (Pcrit * Tcrit) * \
            (-0.37464 * (1 - np.sqrt(Tr)) - 0.5 * 1.54226 * omega / np.sqrt(Tr) + 2 * 0.26992 * omega * np.sqrt(Tr))
    dB_dT = 0.07780 * (R * Tcrit) / Pcrit
    dA_dP = -0.45724 * (1 + 1.54226 * omega - 0.26992 * omega ** 2) * (R * Tcrit) ** 2 / (Pcrit ** 2) / Tcrit
    dB_dP = 0.07780 * (R * Tcrit) / (Pcrit ** 2)
    
    return A, B, dA_dT, dB_dT, dA_dP, dB_dP


def find_real_roots(r):
    N = len(r)
    rr=[]
    RR = np.zeros((N))
    for i in range(N):
        if r[i].imag==0.0:        
            if r[i].real > 0.0:
               RR[i] = r[i]
    for i in range(N):
        if RR[i] > 0.0:
           r_positive = RR[i]
           rr.append(r_positive)
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
        #print ("More than one roots exist for PR equations for H2-N2 Mixtures")
        RR=[]
        RR.append(real_roots[0])
        #print ("RR=", RR)
        phi =  fugacity(RR, A, B)
        real_roots = RR
#        XX  = np.where(real_roots > B )[0] # index where roots are > B
#        real_RR=[]
#        for i in range(len(XX)):
#         #   print ("real_roots[XX]", real_roots[XX[i]])
#            ROOt = real_roots[XX[i]]
#            real_RR.append(ROOt)
#        #print("ROOt", real_RR)
#        phi = fugacity(real_RR, A, B) # find fugacity of those roots
#        print ("phi=",phi)
#        if len(phi) > 1:
#            XXX = np.where(phi == min(phi))[0]
#         #   print ("min fug", XXX[0])
#            real_roots = real_RR[XXX[0]]
#          #  print ("real_root min fug",real_roots)
#            phi = phi[XXX[0]]
    elif len(real_roots)==1:
        #print ("Just one root exist for PR equation for H2-N2 Mixtures")
        phi = fugacity(real_roots, A, B)
    return real_roots, phi



def PR_H2_N2_Main(T, P, x_H2):
    
    #Calculate Coeffs for Peng-Robinson Polynomial
    A, B, alpha, Tr, m = PR_H2_N2_Coeff(T, P, x_H2 )
    
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
    #print ("real_roots are",real_roots)
    ZZ , phi_ZZ = Sort_roots(real_roots, A, B)
    #print("After roots sorting ZZ=", ZZ)
    return ZZ, phi_ZZ


def Cal_Therm(P, T, xH2):
    Tcrit = np.array([33.2, 126.2]) #in kelvin
    Pcrit = np.array([1.297*10**6, 3.394*10**6]) # in Pa
    Tref = xH2*Tcrit[0]+(1-xH2)*Tcrit[1]
    Pref = xH2*Pcrit[0]+(1-xH2)*Pcrit[1]
    xN2 = 1-xH2
    R = 8.3145 #2.016
    Cp_Coeff_N2 = np.array([-2.07E-09/6, 2.19E-06/5, -9.09E-04/4, 1.86E-01/3, -1.89E+01/2, 1.81E+03])
    Cp_Coeff_H2 = np.array([-2.64032606e-15/9,  4.34339102e-12/8, -2.94216657e-09/7, \
                            1.04846126e-06/6, -2.04658359e-04/5,  1.94881968e-02/4, \
                            -3.77290489e-01/3, -3.89881763e+01/2, 1.19295272e+04])
    Cp_Coeff_O2 = np.array([-7.74E-08/6, 7.91E-05/5, -3.14E-02/4, 6.06E+00/3, -5.66E+02/2, 21311.8229])
    H_ideal = (Cp_Coeff_N2[0]*(T**6- Tref**6) + Cp_Coeff_N2[1]*(T**5- Tref**5) + Cp_Coeff_N2[2]*(T**4- Tref**4) + \
              Cp_Coeff_N2[3]*(T**3- Tref**3) + Cp_Coeff_N2[4]*(T**2- Tref**2) + Cp_Coeff_N2[5]*(T- Tref))*xN2  + \
              (Cp_Coeff_H2[0]*(T**9- Tref**9) + Cp_Coeff_H2[1]*(T**8- Tref**8) + Cp_Coeff_H2[2]*(T**7- Tref**7) + \
               Cp_Coeff_H2[3]*(T**6- Tref**6) + Cp_Coeff_H2[4]*(T**5- Tref**5) + Cp_Coeff_H2[5]*(T**4- Tref**4) + \
               Cp_Coeff_H2[6]*(T**3- Tref**3) + Cp_Coeff_H2[7]*(T**2- Tref**2) + Cp_Coeff_H2[8]*(T- Tref))*xH2 
    # Departure from ideal gas behaviour
    Z, phi = PR_H2_N2_Main(T, P, xH2)
    print ("Z for mixtures=", Z)
    A, B, alpha, Tr, m = PR_H2_N2_Coeff(T, P, xH2)
    H_dep = Cal_H_dep(Z, A, B, T, alpha, Tr, R, m)
    print ("Hdep=", H_dep)
    print ("Hideal=", H_ideal)
    Z_ref, phi_ref = PR_H2_N2_Main(Tref, Pref, xH2)
    H_ref = Cal_H_dep(Z_ref, A, B, Tref, alpha, Tr, R, m)
    print ("Href=", H_ref)
    H_dep = xH2*H_dep[0] + xN2*H_dep[1]
    H_ref = xH2*H_ref[0] + xN2*H_ref[1]
    print ("H_dep final=", H_dep)
    print ("H_ref final=", H_ref)
    H = H_dep + H_ideal - H_ref
    print ("H_final=", H)
    return H, H_ref

def Cal_Therm_new(P, T, xH2):
    Tcrit = np.array([33.2, 126.2]) #in kelvin
    Pcrit = np.array([1.297*10**6, 3.394*10**6]) # in Pa
    Tref = xH2*Tcrit[0]+(1-xH2)*Tcrit[1]
    Pref = xH2*Pcrit[0]+(1-xH2)*Pcrit[1]
    xN2 = 1-xH2
    R = 8.3145 #2.016
    Cp_Coeff_N2 = np.array([-2.07E-09/6, 2.19E-06/5, -9.09E-04/4, 1.86E-01/3, -1.89E+01/2, 1.81E+03])
    Cp_Coeff_H2 = np.array([-2.64032606e-15/9,  4.34339102e-12/8, -2.94216657e-09/7, \
                            1.04846126e-06/6, -2.04658359e-04/5,  1.94881968e-02/4, \
                            -3.77290489e-01/3, -3.89881763e+01/2, 1.19295272e+04])
    Cp_Coeff_O2 = np.array([-7.74E-08/6, 7.91E-05/5, -3.14E-02/4, 6.06E+00/3, -5.66E+02/2, 21311.8229])
    H_ideal = (Cp_Coeff_N2[0]*(T**6- Tref**6) + Cp_Coeff_N2[1]*(T**5- Tref**5) + Cp_Coeff_N2[2]*(T**4- Tref**4) + \
              Cp_Coeff_N2[3]*(T**3- Tref**3) + Cp_Coeff_N2[4]*(T**2- Tref**2) + Cp_Coeff_N2[5]*(T- Tref))*xN2  + \
              (Cp_Coeff_H2[0]*(T**9- Tref**9) + Cp_Coeff_H2[1]*(T**8- Tref**8) + Cp_Coeff_H2[2]*(T**7- Tref**7) + \
               Cp_Coeff_H2[3]*(T**6- Tref**6) + Cp_Coeff_H2[4]*(T**5- Tref**5) + Cp_Coeff_H2[5]*(T**4- Tref**4) + \
               Cp_Coeff_H2[6]*(T**3- Tref**3) + Cp_Coeff_H2[7]*(T**2- Tref**2) + Cp_Coeff_H2[8]*(T- Tref))*xH2 
    # Departure from ideal gas behaviour
    Z, phi = PR_H2_N2_Main(T, P, xH2)
    print ("Z for mixtures=", Z)
    A, B, alpha, Tr, m = PR_H2_N2_Coeff(T, P, xH2)
    H_dep = Cal_H_dep_new(Z, A, B, T, alpha, Tr, R, m, xH2, P)
    print ("Hdep=", H_dep)
    print ("Hideal=", H_ideal)
    Z_ref, phi_ref = PR_H2_N2_Main(Tref, Pref, xH2)
    H_ref = Cal_H_dep_new(Z_ref, A, B, Tref, alpha, Tr, R, m , xH2, P)
    print ("Href=", H_ref)
    #H_dep = xH2*H_dep[0] + xN2*H_dep[1]
    #H_ref = xH2*H_ref[0] + xN2*H_ref[1]
    print ("H_dep final=", H_dep)
    print ("H_ref final=", H_ref)
    H = H_dep + H_ideal - H_ref
    print ("H_final=", H)
    return H_dep, H_ref

def Cal_H_dep(Z, A, B, T, alpha, Tr, R, m):
    sqrt2 = np.sqrt(2)
    print ("m=", m)
    print ("alpha=", alpha)
    print ("Tr=", Tr)
    try:
        Hdep = R*T*(Z-1.-A/B/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))* \
                    (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
    except:
        Hdep = R*T*(Z[0]-1.-A/B/2/sqrt2*np.log((Z[0]+(1+sqrt2)*B)/(Z[0]+(1-sqrt2)*B))* \
                    (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
        #Hdep = R*T*(Z-1.-A/B/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))* \
        #            (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
    return Hdep


def Cal_H_dep_new(Z, A, B, T, alpha, Tr, R, m, xH2, P):
    x_H2 = xH2
    x_N2 = 1-x_H2
    Tcrit = np.array([33.2, 126.2]) #in kelvin
    Pcrit = np.array([1.297*10**6, 3.394*10**6]) # in Pa
    Tr    = T/Tcrit
    omega = np.array([-0.219,  0.039])
    R = 8.3145
    m =  0.37464 + 1.54226*omega - 0.26992*omega**2
    alpha = (1+m*(1-np.sqrt(Tr)))**2
    a     = 0.45724*alpha*((R*Tcrit)**2)/Pcrit
    b     = 0.07780*(R*Tcrit)/Pcrit
    k_12  = -0.0878 #from ref 101
    k_11  = 0.0
    k_22  = 0.0
    a11   = a[0]
    a22   = a[1]
    a12   = (1.-k_12)*(a[0]*a[1])**0.5
    b11   = b[0]
    b22   = b[1]
    b12   = (b[0]+b[1])/2.0 # From Clapeyron.jl code
    a_m   = x_H2*x_H2*a11 + x_N2*x_N2*a22 + 2.0*x_H2*x_N2*a12
    b_m   = x_H2*x_H2*b11 + x_N2*x_N2*b22 + x_H2*x_N2*b12
    A     = (a_m*P)/(R*T)**2
    B     = (b_m*P)/(R*T)
    sqrt2 = np.sqrt(2)
    print ("a_m=", a_m)
    print ("alpha=", alpha)
    print ("Z=", Z[0])
    Hdep = R*T*(Z[0]-1.-1.0/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))*\
                -(x_H2*x_H2*a[0]+x_N2*x_N2*a[1]+2.0*x_H2*x_N2*(1-k_12)*(a[0]*a[1])**0.5)-\
                    ((x_H2*x_H2*(a[0]/(Tcrit[0]*Tr[0]**0.5))*(alpha[0])**0.5)+\
                    (x_N2*x_N2*(a[1]/(Tcrit[1]*Tr[1]**0.5))*(alpha[1])**0.5)+\
                    (x_H2*x_N2*(1-k_12)*\
                    ((a[0]**0.5/a[1]**0.5)*a[1]*omega[1]*(alpha[1]**0.5)/Tcrit[1]*(Tr[1]**0.5)+
                     (a[1]**0.5/a[0]**0.5)*a[0]*omega[0]*(alpha[0]**0.5)/Tcrit[0]*(Tr[0]**0.5)))\
                            ))
    #try:
    #    Hdep = R*T*(Z-1.-A/B/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))* \
    #                (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
    #except:
    #    Hdep = R*T*(Z[0]-1.-A/B/2/sqrt2*np.log((Z[0]+(1+sqrt2)*B)/(Z[0]+(1-sqrt2)*B))* \
    #                (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
        #Hdep = R*T*(Z-1.-A/B/2/sqrt2*np.log((Z+(1+sqrt2)*B)/(Z+(1-sqrt2)*B))* \
        #            (1+m*np.sqrt(Tr)/np.sqrt(alpha)))
    return Hdep





def Cal_Cp_mixture(T, xH2):
    Cp_Coeff_N2=np.array([1.50848431e-15, -2.46381879e-12,  1.73643875e-09, -6.89801002e-07, \
                          1.69042289e-04, -2.62077000e-02,  2.51795719e+00, -1.37898047e+02, \
                          4.37824054e+03])
    Cp_Coeff_H2 = np.array([-2.64032606e-15,  4.34339102e-12, -2.94216657e-09, \
                            1.04846126e-06, -2.04658359e-04,  1.94881968e-02, \
                            -3.77290489e-01, -3.89881763e+01, 1.19295272e+04])
    Cp_H2 = Cp_Coeff_H2[0]*(T**8) + Cp_Coeff_H2[1]*(T**7) + Cp_Coeff_H2[2]*(T**6) + \
                Cp_Coeff_H2[3]*(T**5) + Cp_Coeff_H2[4]*(T**4) + Cp_Coeff_H2[5]*(T**3) + \
                Cp_Coeff_H2[6]*(T**2) + Cp_Coeff_H2[7]*(T**1) + Cp_Coeff_H2[8]
    Cp_N2 = (Cp_Coeff_N2[0]*(T**8) + Cp_Coeff_N2[1]*(T**7) + Cp_Coeff_N2[2]*(T**6) + \
             Cp_Coeff_N2[3]*(T**5) + Cp_Coeff_N2[4]*(T**4) + Cp_Coeff_N2[5]*(T**3) + \
             Cp_Coeff_N2[6]*(T**2) + Cp_Coeff_N2[7]*(T**1) + Cp_Coeff_N2[8])
    Cp_mix = Cp_H2*xH2 + (1-xH2)*Cp_N2
    return Cp_mix

import numpy as np

def PR_H2_N2_Derivatives(T, P, x_H2):
    x_N2 = 1 - x_H2
    Tcrit = np.array([33.2, 126.2])  # in Kelvin
    Pcrit = np.array([1.297e6, 3.394e6])  # in Pa
    Tr = T / Tcrit
    omega = np.array([-0.219, 0.039])
    R = 8.3145
    
    # Calculate derivatives of A and B with respect to temperature and pressure
    dA_dT = 0.45724 * (1 + 1.54226 * omega - 0.26992 * omega ** 2) * (R ** 2 * Tcrit ** 2) / (Pcrit * Tcrit) * \
            (-0.37464 * (1 - np.sqrt(Tr)) - 0.5 * 1.54226 * omega / np.sqrt(Tr) + 2 * 0.26992 * omega * np.sqrt(Tr))
    dB_dT = 0.07780 * (R * Tcrit) / Pcrit
    dA_dP = -0.45724 * (1 + 1.54226 * omega - 0.26992 * omega ** 2) * (R * Tcrit) ** 2 / (Pcrit ** 2) / Tcrit
    dB_dP = 0.07780 * (R * Tcrit) / (Pcrit ** 2)
    
    return dA_dT, dB_dT, dA_dP, dB_dP

# Example usage
#T = 300  # Temperature in K
#P = 1e5  # Pressure in Pa
#x_H2 = 0.5  # Mole fraction of H2

#dA_dT, dB_dT, dA_dP, dB_dP = PR_H2_N2_Derivatives(T, P, x_H2)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def find_rho_PT_Mixture(p,T, x_H2): #done
    mol_mass = ((x_H2*2.016 + (1.0-x_H2)*28.010))*10**-3
    R = 8.3145
    Z, phi = PR_H2_N2_Main(T, p, x_H2)
    try:
        rho_val = (p/(R*T))*(mol_mass/Z) 
    except:
        rho_val = (p/(R*T))*(mol_mass/Z[0]) 
    return rho_val
    
def find_H_PT(p,T, x_H2): #done
    H, Href = Cal_Therm(p, T, x_H2)
    return H

def find_H_PT_new(p,T, x_H2): #done
    H, Href = Cal_Therm_new(p, T, x_H2)
    return H[0]

def find_T_rhoP_Mixture(rho, p, xH2): #done
    ## Initial guess for P from ideal gas equation
    R = 8.314462618 #8.3145
    tol = 10.0
    MW_H2 = 2.01588*10**-3
    MW_N2 = 0.02896546 #28.010*10**-3
    Y_ = xH2*MW_H2/(xH2*MW_H2+(1-xH2)*MW_N2)
    mol_mass = MW_N2*MW_H2/(Y_*(MW_N2-MW_H2) + MW_H2)
    #mol_mass = 1/(xH2/MW_H2 + (1-xH2)/MW_N2)
    T = (p*mol_mass)/(rho*R)
    #print ("T=", T)
    #T_cl0 = ambient.P*MW_cl0/(const.R*self.initial_node.rho_cl)
    T1 = T - (tol/100.)*T
    T2 = T + (tol/100.)*T
    #print ("T1=", T1)
    #print ("T2=", T2)
    #print ("xH2=", xH2)
    T_final = regulaFalsi_T(T1, T2, p, rho, xH2)
    return T_final

def derivative_T_rhoP_Mixture(rho, p, xH2):
    # Calculate T using existing function
    MW_H2 = 2.01588*10**-3
    MW_N2 = 0.02896546
    P = 101325.0
    R = 8.314462618 #8.3145
    Y_ = xH2*MW_H2/(xH2*MW_H2+(1-xH2)*MW_N2)
    mol = MW_N2*MW_H2/(Y_*(MW_N2-MW_H2) + MW_H2)
    T = find_T_rhoP_Mixture(rho, p, xH2)
    a, b = PR_H2_N2_Coeff_ab(T, P, xH2)
    M = mol
    dT_drho = 2*a*(mol/R)*((mol**2 - (b*rho)**2)/(mol**2 + 2*mol*b*rho - (b*rho)**2)**2)\
              - ((mol*T)/rho)*(1./(mol-b*rho)) 
    return dT_drho

def derivative_T_Y_Mixture(rho, p, xH2):
    # Calculate T using existing function
    MW_H2 = 2.01588*10**-3
    MW_N2 = 0.02896546 
    R = 8.314462618 #8.3145
    P = 101325.0
    Y_ = xH2*MW_H2/(xH2*MW_H2+(1-xH2)*MW_N2)
    mol = MW_N2*MW_H2/(Y_*(MW_N2-MW_H2) + MW_H2)
    T = find_T_rhoP_Mixture(rho, p, xH2)
    dmoldY = ((MW_H2-MW_N2)*(MW_H2*MW_N2))/(MW_H2 + Y_*(MW_N2-MW_H2))**2
    am, bm = PR_H2_N2_Coeff_ab(T, P, xH2)
    a, b = PR_H2_N2_Coeff_ab_ind(T, P, xH2)
    F_N2 = 0.767
    F_O2 = 0.237
    Ka_New =   [[0.0, 0.0436, 0.433],
                [0.0436, 0.0, -0.0821],
                [0.433, -0.0821, 0.0]]
    Kb_New =   [[0.0, -0.1611, 0.0382],
                [-0.1611, 0.0, -0.2633],
                [0.0382, -0.2633, 0.0]]
    dadY = (2*xH2*a[0]*(1-Ka_New[0][0]) - 2*F_N2*(1-xH2)*a[1]*(1-Ka_New[1][1]) - 2*a[2]*(F_O2**2)*(1-xH2)*(1-Ka_New[2][2])\
            + 2*(1-Ka_New[0][1])*F_N2*(1-2*xH2)*(a[0]*a[1])**0.5 + 2*(1-Ka_New[0][2])*F_O2*(1-2*xH2)*(a[0]*a[2])**0.5 \
            - 4*(1-Ka_New[1][2])*F_N2*F_O2*(1-xH2)*(a[1]*a[2])**0.5)*((MW_H2*MW_N2)/(MW_H2 + Y_*(MW_N2-MW_H2))**2)
    dbdY = (2*xH2*b[0]*(1-Kb_New[0][0]) - 2*F_N2*(1-xH2)*b[1]*(1-Kb_New[1][1]) - 2*F_O2*(1-xH2)*(1-Kb_New[2][2])*b[2] \
            + (1-Kb_New[0][1])*F_N2*(1-2*xH2)*(b[0]+b[1]) + (1-Kb_New[0][2])*F_O2*(1-2*xH2)*(b[0]+b[2])\
                - (1-Kb_New[1][2])*F_N2*F_O2*(1-xH2)*(b[1]+b[2]))*((MW_H2*MW_N2)/(MW_H2 + Y_*(MW_N2-MW_H2))**2)
            
    #        - 2*(1-xH2)*b22 + b12*(1-2*xH2))*((MW_H2*MW_N2)/(MW_H2 + Y_*(MW_N2-MW_H2))**2)
    dTdY = (T/(mol-bm*rho))*(dmoldY - rho*dbdY)\
        + ((rho/R)*((mol-bm*rho)/(mol**2 + 2*mol*bm*rho - (bm*rho)**2)))*(dadY - ((2*am)/(mol**2 + 2*mol*bm*rho - (bm*rho)**2))\
                                                                       *((mol+bm*rho)*dmoldY + (mol-bm*rho)*rho*dbdY))
      #  ((a*rho**2)/((mol**2 + 2*mol*b*rho - (b*rho)**2)**2))*((mol+b*rho)*2*dmoldY + (rho*mol-b*rho**2)*2*dbdY)\
      #  - (dadY*rho**2)/((mol**2 + 2*mol*b*rho - (b*rho)**2))
    return dTdY

def derivative_T_rhoP_Mixture_H(rho, p, xH2):
    # Calculate T using existing function
    T = find_T_rhoP_Mixture(rho, p, xH2)
    
    # Perturb rho slightly for finite difference calculation
    delta_rho = 1e-8 * rho  # Smaller perturbation for higher accuracy
    rho_plus = rho + delta_rho
    rho_minus = rho - delta_rho
    
    # Calculate T for rho+ and rho- using existing function
    T_plus = find_T_rhoP_Mixture(rho_plus, p, xH2)
    T_minus = find_T_rhoP_Mixture(rho_minus, p, xH2)
    
    # Calculate derivative using central difference approximation
    dT_drho = (-T_plus + 8 * T - 8 * T_minus + T_plus) / (12 * delta_rho)
    
    return dT_drho


def dT_drho_PR_H2_N2(T, P, x_H2):
    R = 8.314462618 #8.3145
    # Calculate derivatives of A and B with respect to temperature and pressure
    #dA_dT, dB_dT, dA_dP, dB_dP = PR_H2_N2_Derivatives(T, P, x_H2)
    A, B, dA_dT, dB_dT, dA_dP, dB_dP = PR_H2_N2_Coeff_deriv(T, P, x_H2)
    #dA_dT = dA_dT
    #dB_dT = dB_dT
    #dA_dP = dA_dP
    #dB_dP = dB_dP
    
    # Calculate derivatives of Z with respect to temperature and pressure
    dZ_dT = -(dA_dT - 2*P*dA_dP)/(R*T**2)
    dZ_dP = (B - 2*P*dB_dP)/(R*T)
    Z, phi = PR_H2_N2_Main(T, P, x_H2)
    # Calculate derivatives of T with respect to density using chain rule
    #dT_drho = dT_dZ*dZ_dT + dT_dP*dP_dZ*dZ_dT + dT_dP*dP_dZ*dZ_dP
    try:
        #rho_val = (p/(R*T))*(mol_mass/Z) 
        dT_drho = -T * dZ_dT / (Z * (Z + B))
    except:
        #rho_val = (p/(R*T))*(mol_mass/Z[0])
        dT_drho = -T * dZ_dT / (Z[0] * (Z[0] + B))
    ##dT_drho = -T * dZ_dT / (Z * (Z + B))
    return dT_drho


def find_Cp_H2_N2_Mixture(T, xH2):
    Cp_mix = Cal_Cp_mixture(T, xH2)
    return Cp_mix


def find_T_PH_Mixture(p, H, xH2): #done
    ## Initial guess for P from ideal gas equation
    R = 8.3145
    cp_N2 = 4.37824054e+03
    cp_H2 = 1.19295272e+04
    Tcrit = np.array([33.2, 126.2]) #in kelvin
    Pcrit = np.array([1.297*10**6, 3.394*10**6]) # in Pa
    Tref = xH2*Tcrit[0]+(1-xH2)*Tcrit[1]
    Pref = xH2*Pcrit[0]+(1-xH2)*Pcrit[1]
    #Tref =  33.2*xH2 + 126.2*(1-xH2)
    cp = xH2*cp_H2 + (1-xH2)*cp_N2
    T = H/cp + Tref
    #print ("Assumed values of T=", T)
    tol = 25
    T1 = T - (tol/100.)*T
    T2 = T + (tol/100.)*T
    #print ("T_1=", T1)
    #print ("T_2=", T2)
    T_final = regulaFalsi_H(T1, T2, p, H, xH2)
    print ("T_final=", T_final)
    return T_final



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def func_TPrho(T, p, rho, xH2):
    ## for P
    R = 8.3145
    MW_H2 = 2.01588*10**-3
    MW_N2 = 0.02896546 #28.010*10**-3
    Y_ = xH2*MW_H2/(xH2*MW_H2+(1-xH2)*MW_N2)
    mol_mass = MW_N2*MW_H2/(Y_*(MW_N2-MW_H2) + MW_H2)
    #mol_mass = 1/(xH2/MW_H2 + (1-xH2)/MW_N2)
    Z, phi = PR_H2_N2_Main(T, p, xH2)
    #print ("ZZZ", Z)
    try:
        X = p/Z - rho*R*T/mol_mass
    except:
        X = p/Z[0] - rho*R*T/mol_mass
    return X


def func_H_S(T, P, H, xH2):
    ## for P
    R = 8.3145
    mol_mass = 2.01588*10**-3
    A, B, alpha, Tr, m = PR_H2_N2_Coeff(T, P, xH2)
    Z, phi = PR_H2_N2_Main(T, P, xH2)
    H_g, H_ref= Cal_Therm(P, T, xH2)
    X = H - H_g
    return X


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def regulaFalsi_H( a , b, p, h, xH2):
    MAX_ITER = 1000000
    #print ("func_H_S(a, p, h, xH2)=", func_H_S(a, p, h, xH2))
    #print ("func_H_S(b, p, h, xH2)=", func_H_S(b, p, h, xH2))
    if func_H_S(a, p, h, xH2) * func_H_S(b, p, h, xH2) >= 0:
        print("You have not assumed right a and b  in pure H2 regulaFalsi_S")
        return -1
     
    c = a # Initialize result
     
    for i in range(MAX_ITER):
        # Find the point that touches x axis
        c = (a * func_H_S(b, p, h, xH2) - b * func_H_S(a, p, h, xH2))/ (func_H_S(b, p, h, xH2) - func_H_S(a, p, h, xH2))
         
        # Check if the above found point is root
        if func_H_S(c, p, h, xH2) == 0:
            break
         
        # Decide the side to repeat the steps
        elif func_H_S(c, p, h, xH2) * func_H_S(a, p, h, xH2) < 0:
            b = c
        else:
            a = c
        return c
  
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