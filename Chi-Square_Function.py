import numpy as np
import cmath
from cmath import sqrt

# Quarck masses
mu = 1.23 
md = 2.67 
ms = 53.16 
mc = 620 
mb = 2839 
mt = 168260 

# Theoretical values of the matrix elements
vud_th = 0.97401 
vus_th = 0.22650 
vub_th = 0.00361 
vcd_th = 0.22636 
vcs_th = 0.97320 
vcb_th = 0.04053 
vtd_th = 0.00854 
vts_th = 0.03978 
vtb_th = 0.999172 
jarslkog_th = 0.00003  

# Uncertainties 
sigma_vud = 0.00011 
sigma_vus = 0.00048 
sigma_vub = 0.00009 
sigma_vcd = 0.00048 
sigma_vcs = 0.00011 
sigma_vcb = 0.00061 
sigma_vtd = 0.00016 
sigma_vts = 0.00060 
sigma_vtb = 0.000024 
sigma_jars = 0.0000009 

# Elements of the matrix O_u that diagonalize M_u
eta_u = 1
eta_d = 1

def O11_u(A, m1, m2, m3):

    o11 = (sqrt(((m3 * (eta_u*m2)) * (A - (-eta_u*m1))) / (A * ((eta_u*m2) - (-eta_u*m1)) * (m3 - (-eta_u*m1)))))
    return o11

def O22_u(A, m1, m2, m3):
    o22 = (sqrt(((eta_u*m2) * (A - (eta_u*m2))) / (((eta_u*m2) - (-eta_u*m1)) * (m3 - (eta_u*m2)))))
    return o22

def O33_u(A, m1, m2, m3):
    o33 = (sqrt((m3 * (A - (-eta_u*m1)) * (A - (eta_u*m2))) / (A * (m3 - (-eta_u*m1)) * (m3 - (eta_u*m2)))))
    return o33

# Off-diagonal elements 

def O12_u(A, m1, m2, m3):
    o12 = (eta_u) * (sqrt(complex((((-eta_u*m1) * m3) * ((eta_u*m2) - A)) / (A * ((eta_u*m2) - (-eta_u*m1)) * (m3 - (eta_u*m2))))))
    return o12

def O13_u(A, m1, m2, m3):
    o13 = (sqrt(complex((((-eta_u*m1) * (eta_u*m2)) * (A - m3)) / (A * (m3 - (-eta_u*m1)) * (m3 - (eta_u*m2))))))
    return o13

def O21_u(A, m1, m2, m3):
    o21 = (-eta_u) * (sqrt(complex(((-eta_u*m1) * ((-eta_u*m1) - A)) / (((eta_u*m2) - (-eta_u*m1)) * (m3 - (-eta_u*m1))))))
    return o21

def O23_u(A, m1, m2, m3):
    o23 = (sqrt(complex((m3 * (m3 - A)) / ((m3 - (-eta_u*m1)) * (m3 - (eta_u*m2))))))
    return o23

def O31_u(A, m1, m2, m3):
    o31 = (eta_u) * (sqrt(complex(((-eta_u*m1) * (A - (eta_u*m2)) * (A - m3)) / (A * ((eta_u*m2) - (-eta_u*m1)) * (m3 - (-eta_u*m1))))))
    return o31

def O32_u(A, m1, m2, m3):
    o32 = (-1) * (sqrt(complex(((eta_u*m2) * (A - (-eta_u*m1)) * (m3 - A)) / (A * ((eta_u*m2) - (-eta_u*m1)) * (m3 - (eta_u*m2))))))
    return o32

# Elements of the matrix O_d that diagonalize M_d

def O11_d(A, m1, m2, m3):

    o11 = (sqrt(((m3 * (eta_d*m2)) * (A - (-eta_d*m1))) / (A * ((eta_d*m2) - (-eta_d*m1)) * (m3 - (-eta_d*m1)))))
    return o11

def O22_d(A, m1, m2, m3):
    o22 = (sqrt(((eta_d*m2) * (A - (eta_d*m2))) / (((eta_d*m2) - (-eta_d*m1)) * (m3 - (eta_d*m2)))))
    return o22

def O33_d(A, m1, m2, m3):
    o33 = (sqrt((m3 * (A - (-eta_d*m1)) * (A - (eta_d*m2))) / (A * (m3 - (-eta_d*m1)) * (m3 - (eta_d*m2)))))
    return o33

# Off-diagonal elements 

def O12_d(A, m1, m2, m3):
    o12 = (eta_d) * (sqrt(complex((((-eta_d*m1) * m3) * ((eta_d*m2) - A)) / (A * ((eta_d*m2) - (-eta_d*m1)) * (m3 - (eta_d*m2))))))
    return o12

def O13_d(A, m1, m2, m3):
    o13 = (sqrt(complex((((-eta_d*m1) * (eta_d*m2)) * (A - m3)) / (A * (m3 - (-eta_d*m1)) * (m3 - (eta_d*m2))))))
    return o13

def O21_d(A, m1, m2, m3):
    o21 = (-eta_d) * (sqrt(complex(((-eta_d*m1) * ((-eta_d*m1) - A)) / (((eta_d*m2) - (-eta_d*m1)) * (m3 - (-eta_d*m1))))))
    return o21

def O23_d(A, m1, m2, m3):
    o23 = (sqrt(complex((m3 * (m3 - A)) / ((m3 - (-eta_d*m1)) * (m3 - (eta_d*m2))))))
    return o23

def O31_d(A, m1, m2, m3):
    o31 = (eta_d) * (sqrt(complex(((-eta_d*m1) * (A - (eta_d*m2)) * (A - m3)) / (A * ((eta_d*m2) - (-eta_d*m1)) * (m3 - (-eta_d*m1))))))
    return o31

def O32_d(A, m1, m2, m3):
    o32 = (-1) * (sqrt(complex(((eta_d*m2) * (A - (-eta_d*m1)) * (m3 - A)) / (A * ((eta_d*m2) - (-eta_d*m1)) * (m3 - (eta_d*m2))))))
    return o32

# Theoretical elements of the Vckm matrix 

def CKM11(Au, Ad, phi1, phi2):
    ckm11 = np.clongdouble((O11_u(Au, mu, mc, mt)*O11_d(Ad, md, ms, mb)) + (O21_u(Au, mu, mc, mt)*O21_d(Ad, md, ms, mb)*cmath.exp(complex(0,phi1))) + (O31_u(Au, mu, mc, mt)*O31_d(Ad, md, ms, mb)*cmath.exp(complex(0,(phi1+phi2)))))
    return ckm11

def CKM12(Au, Ad, phi1, phi2):
    ckm12 = np.clongdouble((O11_u(Au, mu, mc, mt)*O12_d(Ad, md, ms, mb)) + (O21_u(Au, mu, mc, mt)*O22_d(Ad, md, ms, mb)*cmath.exp(complex(0,phi1))) + (O31_u(Au, mu, mc, mt)*O32_d(Ad, md, ms, mb)*cmath.exp(complex(0,(phi1+phi2)))))
    return ckm12

def CKM13(Au, Ad, phi1, phi2):
    ckm13 = np.clongdouble((O11_u(Au, mu, mc, mt)*O13_d(Ad, md, ms, mb)) + (O21_u(Au, mu, mc, mt)*O23_d(Ad, md, ms, mb)*cmath.exp(complex(0,phi1))) + (O31_u(Au, mu, mc, mt)*O33_d(Ad, md, ms, mb)*cmath.exp(complex(0,(phi1+phi2)))))
    return ckm13 

def CKM21(Au, Ad, phi1, phi2):
    ckm21=np.clongdouble((O12_u(Au, mu, mc, mt)*O11_d(Ad, md, ms, mb)) + (O22_u(Au, mu, mc, mt)*O21_d(Ad, md, ms, mb)*cmath.exp(complex(0,phi1))) + (O32_u(Au, mu, mc, mt)*O31_d(Ad, md, ms, mb)*cmath.exp(complex(0,(phi1+phi2)))))
    return ckm21

def CKM22(Au, Ad, phi1, phi2):
    ckm22=np.clongdouble((O12_u(Au, mu, mc, mt)*O12_d(Ad, md, ms, mb)) + (O22_u(Au, mu, mc, mt)*O22_d(Ad, md, ms, mb)*cmath.exp(complex(0,phi1))) + (O32_u(Au, mu, mc, mt)*O32_d(Ad, md, ms, mb)*cmath.exp(complex(0,(phi1+phi2)))))
    return ckm22

def CKM23(Au, Ad, phi1, phi2):
    ckm23=np.clongdouble((O12_u(Au, mu, mc, mt)*O13_d(Ad, md, ms, mb)) + (O22_u(Au, mu, mc, mt)*O23_d(Ad, md, ms, mb)*cmath.exp(complex(0,phi1))) + (O32_u(Au, mu, mc, mt)*O33_d(Ad, md, ms, mb)*cmath.exp(complex(0,(phi1+phi2)))))
    return ckm23

def CKM31(Au, Ad, phi1, phi2):
    ckm31=np.clongdouble((O13_u(Au, mu, mc, mt)*O11_d(Ad, md, ms, mb)) + (O23_u(Au, mu, mc, mt)*O21_d(Ad, md, ms, mb)*cmath.exp(complex(0,phi1))) + (O33_u(Au, mu, mc, mt)*O31_d(Ad, md, ms, mb)*cmath.exp(complex(0,(phi1+phi2)))))
    return ckm31
    
def CKM32(Au, Ad, phi1, phi2):
    ckm32=np.clongdouble((O13_u(Au, mu, mc, mt)*O12_d(Ad, md, ms, mb)) + (O23_u(Au, mu, mc, mt)*O22_d(Ad, md, ms, mb)*cmath.exp(complex(0,phi1))) + (O33_u(Au, mu, mc, mt)*O32_d(Ad, md, ms, mb)*cmath.exp(complex(0,(phi1+phi2)))))
    return ckm32

def CKM33(Au, Ad, phi1, phi2):
    ckm33=np.clongdouble((O13_u(Au, mu, mc, mt)*O13_d(Ad, md, ms, mb)) + (O23_u(Au, mu, mc, mt)*O23_d(Ad, md, ms, mb)*cmath.exp(complex(0,phi1))) + (O33_u(Au, mu, mc, mt)*O33_d(Ad, md, ms, mb)*cmath.exp(complex(0,(phi1+phi2)))))
    return ckm33

def JJ(Au, Ad, phi1, phi2):
    jj=np.imag(CKM12(Au, Ad, phi1, phi2)*CKM23(Au, Ad, phi1, phi2)*(np.conj(CKM13(Au, Ad, phi1, phi2)))*(np.conj(CKM22(Au, Ad, phi1, phi2))))
    return jj

# Chi-Square function

def chi2CKM23(Au, Ad, phi1, phi2):
    chi23 = ((vcb_th - abs(CKM23(Au, Ad, phi1, phi2)))**2) / (sigma_vcb**2)
    return chi23

def chi2CKM12(Au, Ad, phi1, phi2):
    chi12 = ((vus_th - abs(CKM12(Au, Ad, phi1, phi2)))**2) / (sigma_vus**2)
    return chi12

def chi2CKM13(Au, Ad, phi1, phi2):
    chi13 = ((vub_th - abs(CKM13(Au, Ad, phi1, phi2)))**2) / (sigma_vub**2)
    return chi13

def chi2JJ( Au,  Ad, phi1, phi2):
    chijj = ((jarslkog_th - JJ(Au, Ad, phi1, phi2))**2) / (sigma_jars**2)
    return chijj

def Chi_Square(Au, Ad, phi1, phi2):
    chi_square = (chi2CKM23(Au, Ad, phi1, phi2) + chi2CKM12(Au, Ad, phi1, phi2) + chi2CKM13(Au, Ad, phi1, phi2) + chi2JJ(Au, Ad, phi1, phi2)) / 4
    return chi_square