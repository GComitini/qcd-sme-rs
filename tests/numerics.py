#!/usr/bin/env python3

import numpy as np


def en(q, m):
    return np.sqrt(q**2+m**2)


def bose(q, m, beta):
    return 1/(np.exp(beta*en(q, m))-1)


def fermi(q, m, beta, mu):
    return 1/(np.exp(beta*(en(q, m)-mu))+1)


def fermi_double(q, m, beta, mu):
    return fermi(q, m, beta, mu)+fermi(q, m, beta, -mu)


def j_m_i(q, m, beta):
    return q**2*bose(q, m, beta)/(en(q, m)*2*np.pi**2)


def d_j_m_i(q, m, beta):
    return -bose(q, m, beta)/(en(q, m)*4*np.pi**2)


def j_m_l_i(q, m, beta):
    return -q**2*bose(q, m, beta)*en(q, m)/(2*np.pi**2)


def j_m_t_i(q, m, beta):
    return q**4*bose(q, m, beta)/(en(q, m)*6*np.pi**2)


def j_0_t_i(q, beta):
    return j_m_t_i(q, 0, beta)


def tlog(q, om, p, en, delta):
    r = om**2+p**2+2j*en*om+delta
    r0 = 2*q*p
    ru = r+r0
    rd = r-r0
    return np.log(ru/rd)


def tlog_1o(q, om, en, delta):
    r = om**2+2j*en*om+delta
    return 4*q/r


def tlog_l(q, om, p, m1, m2):
    en1 = np.sqrt(q*q+m1*m1)
    en2 = np.sqrt(q*q+m2*m2)
    z = om+1j*en2
    delta = m1*m1-m2*m2
    return (z**2+en1**2+p**2*(2*z/om-1))**2*tlog(q, om, p, en2, delta)


def tlog_l_3o(q, om, en, delta):
    r = om**2+2j*en*om+delta
    return 4*q*(1+4j*en/om+4*q**2/(3*r))


def tlog_l2(q, om, p, en, delta):
    return (om**2+2j*om*en+delta+p**2*(2j*en/om+1))**2*tlog(q, om, p, en, delta)


def tlog_l3(q, om, p, en, delta):
    m12 = delta+en**2-q**2
    en1 = np.sqrt(q*q+m12)
    z = om+1j*en
    return (z**2+en1**2+p**2*(2*z/om-1))**2*tlog(q, om, p, en, delta)


def tlog_t(q, om, p, en, delta):
    r = om**2+p**2+2j*en*om+delta
    r0 = 2*q*p
    ru = r+r0
    rd = r-r0
    return ru*rd*tlog(q, om, p, en, delta)


def tlog_t_3o(q, om, en, delta):
    r = om**2+2j*en*om+delta
    return (4*q*(1-8*q**2/(3*r)))


def i_m_m_i(q, om, p, m1, m2, beta):
    delta1 = m2**2-m1**2
    delta2 = m1**2-m2**2
    t1 = bose(q, m1, beta)/en(q, m1)*(tlog(q, om, p, en(q, m1),
                                           delta1)+tlog(q, -om, p, en(q, m1), delta1))/2
    t2 = bose(q, m2, beta)/en(q, m2)*(tlog(q, om, p, en(q, m2),
                                           delta2)+tlog(q, -om, p, en(q, m2), delta2))/2
    return (t1+t2)*q/p/(8*np.pi**2)


def i_m_m_i_zp(q, om, m1, m2, beta):
    delta1 = m2**2-m1**2
    delta2 = m1**2-m2**2
    t1 = bose(q, m1, beta)/en(q, m1)*(tlog_1o(q, om, en(q, m1),
                                              delta1)+tlog_1o(q, -om, en(q, m1), delta1))/2
    t2 = bose(q, m2, beta)/en(q, m2)*(tlog_1o(q, om, en(q, m2),
                                              delta2)+tlog_1o(q, -om, en(q, m2), delta2))/2
    return (t1+t2)*q/(8*np.pi**2)


def i_m_m_l_i(q, om, p, m1, m2, beta):
    om = om + 1E-16
    delta1 = m2**2-m1**2
    delta2 = m1**2-m2**2
    t1 = bose(q, m1, beta)/en(q, m1)*((tlog_l(q, om, p, m2, m1) +
                                       tlog_l(q, -om, p, m2, m1))/2-4*q*p*(om**2+p**2+delta1))
    t2 = bose(q, m2, beta)/en(q, m2)*((tlog_l(q, om, p, m1, m2) +
                                       tlog_l(q, -om, p, m1, m2))/2-4*q*p*(om**2+p**2+delta2))
    return om**2/(om**2+p**2)*(t1+t2)*q/(p**3*32*np.pi**2)


def i_m_m_l_i_zp(q, om, m1, m2, beta):
    delta = m2**2-m1**2
    en1 = en(q, m1)
    en2 = en(q, m2)
    r1p = 1/(om*om+2j*om*en1+delta)
    r1m = 1/(om*om-2j*om*en1+delta)
    r2p = 1/(om*om+2j*om*en2-delta)
    r2m = 1/(om*om-2j*om*en2-delta)
    t1 = bose(q, m1, beta)/en1*(r1p+r1m)/2
    t2 = bose(q, m2, beta)/en2*(r2p+r2m)/2
    return (t1+t2)*q**4/(6*np.pi**2)


def i_m_m_t_i(q, om, p, m1, m2, beta):
    delta1 = m2**2-m1**2
    delta2 = m1**2-m2**2
    en1 = en(q, m1)
    en2 = en(q, m2)
    t1 = bose(q, m1, beta)/en1*((tlog_t(q, om, p, en1, delta1) +
                                 tlog_t(q, -om, p, en1, delta1))/2-4*q*p*(om**2+p**2+delta1))
    t2 = bose(q, m2, beta)/en(q, m2)*((tlog_t(q, om, p, en2, delta2) +
                                       tlog_t(q, -om, p, en2, delta2))/2-4*q*p*(om**2+p**2+delta2))
    return -(t1+t2)*q/(p**3*64*np.pi**2)


def i_m_m_t_i_zp(q, om, m1, m2, beta):
    delta1 = m2**2-m1**2
    delta2 = m1**2-m2**2
    en1 = en(q, m1)
    en2 = en(q, m2)
    r1p = 1/(om*om+2j*om*en1+delta1)
    r1m = 1/(om*om-2j*om*en1+delta1)
    r2p = 1/(om*om+2j*om*en2+delta2)
    r2m = 1/(om*om-2j*om*en2+delta2)
    t1 = bose(q, m1, beta)/en1*(r1p+r1m)/2
    t2 = bose(q, m2, beta)/en2*(r2p+r2m)/2
    return (t1+t2)*q**4/(6*np.pi**2)


def d_j_m_i(q, m, beta):
    return -bose(q, m, beta)/(en(q, m)*4*np.pi**2)


def d2_j_m_i(q, m, beta):
    bos = bose(q, m, beta)
    e = en(q, m)
    bosm = 1/(np.exp(-e*beta)-1)
    return (bos/e**3-beta*bos*bosm/e**2)/(8*np.pi**2)


def d_j_m_l_i(q, m, beta):
    return en(q, m)*bose(q, m, beta)/(4*np.pi**2)


def jnnm(q, m, beta):
    bos = bose(q, m, beta)
    e = en(q, m)
    bosm = 1/(np.exp(-e*beta)-1)
    return bos*bosm/(8*np.pi**2)


def d2_j_m_l_i(q, m, beta):
    return -d_j_m_i(q, m, beta)/2+beta*jnnm(q, m, beta)


def d_j_m_t_i(q, m, beta):
    return -j_m_i(q, m, beta)/2


def d2_j_m_t_i(q, m, beta):
    return -d_j_m_i(q, m, beta)/2


def d_i_m_m_i(q, om, p, m1, m2, beta):
    x1 = -bose(q, m1, beta)/en(q, m1) * \
        (1/((om+1j*en(q, m1))**2+en(q+p, m2)**2) +
         1/((om+1j*en(q, m1))**2+en(q-p, m2)**2) +
         1/((-om+1j*en(q, m1))**2+en(q+p, m2)**2) +
         1/((-om+1j*en(q, m1))**2+en(q-p, m2)**2))/2
    x2 = bose(q, m2, beta)/en(q, m2) * \
        (1/((om+1j*en(q, m2))**2+en(q+p, m1)**2) -
         1/((om+1j*en(q, m2))**2+en(q-p, m1)**2) +
         1/((-om+1j*en(q, m2))**2+en(q+p, m1)**2) -
         1/((-om+1j*en(q, m2))**2+en(q-p, m1)**2))/2
    x3 = -bose(q, m1, beta)/en(q, m1) * \
        (1/((om+1j*en(q, m1))**2+en(q+p, m2)**2) -
         1/((om+1j*en(q, m1))**2+en(q-p, m2)**2) +
         1/((-om+1j*en(q, m1))**2+en(q+p, m2)**2) -
         1/((-om+1j*en(q, m1))**2+en(q-p, m2)**2))/2
    return (x1+(x2+x3)*q/p)/(8*np.pi**2)


def d_i_m_m_i_zp(q, om, m1, m2, beta):
    x1 = -bose(q, m1, beta)/en(q, m1) * \
        (1/((om+1j*en(q, m1))**2+en(q, m2)**2) +
         1/((-om+1j*en(q, m1))**2+en(q, m2)**2))/2
    x2 = bose(q, m2, beta)/en(q, m2) * \
        (1/((om+1j*en(q, m2))**2+en(q, m1)**2)**2 +
         1/((-om+1j*en(q, m2))**2+en(q, m1)**2)**2)/2
    x3 = -bose(q, m1, beta)/en(q, m1) * \
        (1/((om+1j*en(q, m1))**2+en(q, m2)**2)**2 +
         1/((-om+1j*en(q, m1))**2+en(q, m2)**2)**2)/2
    return (x1-2*(x2+x3)*q*q)/(4*np.pi**2)


def d_i_m_m_l_i(q, om, p, m1, m2, beta):
    om = om+1E-16
    bose1 = bose(q, m1, beta)
    en1 = en(q, m1)
    z1p = om+1j*en1
    z1m = -om+1j*en1
    bose2 = bose(q, m2, beta)
    en2 = en(q, m2)
    z2p = om+1j*en2
    z2m = -om+1j*en2
    t1 = -bose1/en1*((1/(z1p**2+(p+q)**2+m2**2)+1/(z1p**2+(p-q)**2+m2**2))*(z1p**2+q**2+m2**2+p**2*(2*z1p/om-1))**2 +
                     (1/(z1m**2+(p+q)**2+m2**2)+1/(z1m**2+(p-q)**2+m2**2))*(z1m**2+q**2+m2**2+p**2*(-2*z1m/om-1))**2)/4
    t2 = bose1/en1*(om**2+p**2+m2**2-m1**2)
    t3 = 2*q**2*(bose1/en1-bose2/en2)
    t4 = q/p*(bose2/en2*((z2p**2+q**2+m1**2+p**2*(2*z2p/om-1))
              * tlog(q, om, p, en2, m1**2-m2**2)+(z2m**2+q**2+m1**2+p**2*(-2*z2m/om-1))
              * tlog(q, -om, p, en2, m1**2-m2**2))
              - bose1/en1*((z1p**2+q**2+m2**2+p**2*(2*z1p/om-1))
              * tlog(q, om, p, en1, m2**2-m1**2)+(z1m**2+q**2+m2**2+p**2*(-2*z1m/om-1))
              * tlog(q, -om, p, en1, m2**2-m1**2)))/2
    t5 = q/p*(bose2/en2*((z2p**2+q**2+m1**2+p**2*(2*z2p/om-1))**2
              * (1/(z2p**2+en(q+p, m1)**2)-1/(z2p**2+en(q-p, m1)**2))+(z2m**2+q**2+m1**2+p**2*(-2*z2m/om-1))**2
              * (1/(z2m**2+en(q+p, m1)**2)-1/(z2m**2+en(q-p, m1)**2)))
              - bose1/en1*((z1p**2+q**2+m2**2+p**2*(2*z1p/om-1))**2
              * (1/(z1p**2+en(q+p, m2)**2)-1/(z1p**2+en(q-p, m2)**2))+(z1m**2+q**2+m2**2+p**2*(-2*z1m/om-1))**2
              * (1/(z1m**2+en(q+p, m2)**2)-1/(z1m**2+en(q-p, m2)**2))))/4
    return om**2/(om**2+p**2)*(t1+t2+t3+t4+t5)/(16*np.pi**2*p**2)


def d_i_m_m_l_i_zp(q, om, m1, m2, beta):
    bose1 = bose(q, m1, beta)
    en1 = en(q, m1)
    bose2 = bose(q, m2, beta)
    en2 = en(q, m2)
    t10 = bose1/en1*(1/(om*om+2j*en1*om+m2*m2-m1*m1) +
                     1/(om*om-2j*en1*om+m2*m2-m1*m1))/2
    t1 = bose1/en1*(1/(om*om+2j*en1*om+m2*m2-m1*m1)**2 +
                    1/(om*om-2j*en1*om+m2*m2-m1*m1)**2)/2
    t2 = bose2/en2*(1/(om*om+2j*en2*om+m1*m1-m2*m2)**2 +
                    1/(om*om-2j*en2*om+m1*m1-m2*m2)**2)/2
    return (-t10/4+(t1-t2)*q**2/6)*q**2/(np.pi**2)


def d_i_m_m_t_i(q, om, p, m1, m2, beta):
    delta = m2**2-m1**2
    bose1 = bose(q, m1, beta)
    en1 = en(q, m1)
    bose2 = bose(q, m2, beta)
    en2 = en(q, m2)
    t1 = -q/p*bose1/en1*(tlog(q, om, p, en1, delta) +
                         tlog(q, -om, p, en1, delta))/(32*np.pi**2)
    t2 = q**2/p**2*(bose2/en2-bose1/en1)/(8*np.pi**2)
    t31 = bose1/en1*(((om+1j*en1)**2+m2**2+q**2+p**2)*tlog(q, om, p, en1, delta) +
                     ((-om+1j*en1)**2+m2**2+q**2+p**2)*tlog(q, -om, p, en1, delta))/2
    t32 = bose2/en2*(((om+1j*en2)**2+m1**2+q**2+p**2)*tlog(q, om, p, en2, -delta) +
                     ((-om+1j*en2)**2+m1**2+q**2+p**2)*tlog(q, -om, p, en2, -delta))/2
    t3 = q/p**3*(t31-t32)/(32*np.pi**2)
    return t1+t2+t3


def ghost_se(q, om, p, m, beta):
    s = om*om+p*p
    m2 = m*m
    a = s+m2
    t1 = -s*a/(2*m2)*i_m_m_i(q, om, p, m, 0, beta)
    t2 = s*s/(2*m2)*i_m_m_i(q, om, p, 0, 0, beta)
    t3 = (2*s-m2)/(4*m2)*(j_m_i(q, m, beta)-j_m_i(q, 0, beta))
    t4 = a*a/4 * d_i_m_m_i(q, om, p, m, 0, beta)
    t5 = -(s-m2)/4*d_j_m_i(q, m, beta)
    return t1+t2+t3+t4+t5


def ghost_se_zp(q, om, m, beta):
    s = om*om
    m2 = m*m
    a = s+m2
    t1 = -s*a/(2*m2)*i_m_m_i_zp(q, om,  m, 0, beta)
    t2 = s*s/(2*m2)*i_m_m_i_zp(q, om,  0, 0, beta)
    t3 = (2*s-m2)/(4*m2)*(j_m_i(q, m, beta)-j_m_i(q, 0, beta))
    t4 = a*a/4 * d_i_m_m_i_zp(q, om, m, 0, beta)
    t5 = -(s-m2)/4*d_j_m_i(q, m, beta)
    return t1+t2+t3+t4+t5


def j_m_l_p_i(q, om, p, m, beta):
    return (p*p*j_m_l_i(q, m, beta)+om*om*j_m_t_i(q, m, beta))/(om*om+p*p)


def d_j_m_l_p_i(q, om, p, m, beta):
    return (p*p*d_j_m_l_i(q, m, beta)+om*om*d_j_m_t_i(q, m, beta))/(om*om+p*p)


def d2_j_m_l_p_i(q, om, p, m, beta):
    return (p*p*d2_j_m_l_i(q, m, beta)+om*om*d2_j_m_t_i(q, m, beta))/(om*om+p*p)


def j_m_l_p_i_zp(q, om, m, beta):
    return j_m_t_i(q, m, beta)


def d_j_m_l_p_i_zp(q, om, m, beta):
    return d_j_m_t_i(q, m, beta)


def d2_j_m_l_p_i_zp(q, om, m, beta):
    return d2_j_m_t_i(q, m, beta)


def gluon_pol_l(q, om, p, m, beta):
    s = om*om+p*p
    s2 = s*s
    m2 = m*m
    m4 = m2*m2

    t1 = (3*s2/(2*m4)-1)*i_m_m_l_i(q, om, p, 0, 0, beta)
    t2 = (4+(3*s2+8*m2*s+4*m4)/(2*m4))*i_m_m_l_i(q, om, p, m, m, beta)
    t3 = -(3*s2+4*m2*s+m4)/m4 * i_m_m_l_i(q, om, p, m, 0, beta)
    t4 = 2*s*(s+2*m2)/m2*i_m_m_i(q, om, p, m, m, beta)
    t5 = -2*s*(s+m2)/m2*i_m_m_i(q, om, p, m, 0, beta)
    t6 = -(2*s+3*m2)/m2*j_m_i(q, m, beta)
    t7 = (2*s+m2)/m2*j_m_i(q, 0, beta)
    t8 = -(8*m2+(s+2*m2)**2/m2)*d_i_m_m_l_i(q, om, p, m, m, beta)
    t9 = (s+m2)**2/m2 * d_i_m_m_l_i(q, om, p, m, 0, beta)
    t10 = -2*s*(s+4*m2)*d_i_m_m_i(q, om, p, m, m, beta)
    t11 = (s+m2)**2*d_i_m_m_i(q, om, p, m, 0, beta)
    t12 = (s+3*m2)*d_j_m_i(q, m, beta)

    t13 = -m4*d2_j_m_i(q, m, beta)
    t14 = (j_m_l_p_i(q, om, p, m, beta)-j_m_l_p_i(q, om, p, 0, beta))/m2
    t15 = -d_j_m_l_p_i(q, om, p, m, beta)
    t16 = m2/2*d2_j_m_l_p_i(q, om, p, m, beta)

    return t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16


def gluon_pol_l_zp(q, om, m, beta):
    s = om*om
    s2 = s*s
    m2 = m*m
    m4 = m2*m2

    t1 = (3*s2/(2*m4)-1)*i_m_m_l_i_zp(q, om, 0, 0, beta)
    t2 = (4+(3*s2+8*m2*s+4*m4)/(2*m4))*i_m_m_l_i_zp(q, om, m, m, beta)
    t3 = -(3*s2+4*m2*s+m4)/m4 * i_m_m_l_i_zp(q, om, m, 0, beta)
    t4 = 2*s*(s+2*m2)/m2*i_m_m_i_zp(q, om, m, m, beta)
    t5 = -2*s*(s+m2)/m2*i_m_m_i_zp(q, om,  m, 0, beta)
    t6 = -(2*s+3*m2)/m2*j_m_i(q, m, beta)
    t7 = (2*s+m2)/m2*j_m_i(q, 0, beta)
    t8 = -(8*m2+(s+2*m2)**2/m2)*d_i_m_m_l_i_zp(q, om, m, m, beta)
    t9 = (s+m2)**2/m2 * d_i_m_m_l_i_zp(q, om, m, 0, beta)
    t10 = -2*s*(s+4*m2)*d_i_m_m_i_zp(q, om, m, m, beta)
    t11 = (s+m2)**2*d_i_m_m_i_zp(q, om, m, 0, beta)
    t12 = (s+3*m2)*d_j_m_i(q, m, beta)

    t13 = -m4*d2_j_m_i(q, m, beta)
    t14 = (j_m_l_p_i_zp(q, om, m, beta)-j_m_l_p_i_zp(q, om, 0, beta))/m2
    t15 = -d_j_m_l_p_i_zp(q, om, m, beta)
    t16 = m2/2*d2_j_m_l_p_i_zp(q, om, m, beta)

    return t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16


def gluon_pol_t(q, om, p, m, beta):
    s = om*om+p*p
    s2 = s*s
    m2 = m*m
    m4 = m2*m2

    t1 = (3*s2/(2*m4)-1)*i_m_m_t_i(q, om, p, 0, 0, beta)
    t2 = (4+(3*s2+8*m2*s+4*m4)/(2*m4))*i_m_m_t_i(q, om, p, m, m, beta)
    t3 = -(3*s2+4*m2*s+m4)/m4 * i_m_m_t_i(q, om, p, m, 0, beta)
    t4 = 2*s*(s+2*m2)/m2*i_m_m_i(q, om, p, m, m, beta)
    t5 = -2*s*(s+m2)/m2*i_m_m_i(q, om, p, m, 0, beta)
    t6 = -(2*s+3*m2)/m2*j_m_i(q, m, beta)
    t7 = (2*s+m2)/m2*j_m_i(q, 0, beta)
    t8 = -(8*m2+(s+2*m2)**2/m2)*d_i_m_m_t_i(q, om, p, m, m, beta)
    t9 = (s+m2)**2/m2 * d_i_m_m_t_i(q, om, p, m, 0, beta)
    t10 = -2*s*(s+4*m2)*d_i_m_m_i(q, om, p, m, m, beta)
    t11 = (s+m2)**2*d_i_m_m_i(q, om, p, m, 0, beta)
    t12 = (s+3*m2)*d_j_m_i(q, m, beta)

    t13 = -m4*d2_j_m_i(q, m, beta)
    t14 = (j_m_t_i(q, m, beta)-j_m_t_i(q, 0, beta))/m2
    t15 = -d_j_m_t_i(q, m, beta)
    t16 = m2/2*d2_j_m_t_i(q, m, beta)

    return t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16


def gluon_pol_t_zp(q, om, m, beta):
    return gluon_pol_l_zp(q, om, m, beta)


def gluon_pol_quark_l(q, om, p, m, beta, mu):
    s = om*om+p*p
    eg = en(q, m)
    t1 = (s-4*eg**2+4j*eg*om)/(8*q*p)*tlog(q, om, p, eg, 0)
    t1opp = (s-4*eg**2-4j*eg*om)/(8*q*p)*tlog(q, -om, p, eg, 0)
    return -(s/(p*p))*fermi_double(q, m, beta, mu)*q**2/eg*(1-t1-t1opp)/(6*np.pi**2)


def gluon_pol_quark_l_zp(q, om, m, beta, mu):
    eg = en(q, m)
    return -2/(9*np.pi**2)*fermi_double(q, m, beta, mu)*q**2/eg*(3*eg**2-q**2)/(om**2+4*eg**2)


Qlog = [
    (0.03, 0., 2.12, 4., 2.),
    (0.047, 0., 2.12, 4., 2.),
    (0.047, 0.18, 2.12, 4., 2.),
    (0.047, 0.18, 1.27, 4., 2.),
    (0.047, 0.18, 1.27, 6.9, 2.),
    (0.047, 0.18, 1.27, 6.9, 3.21),
]

Qlogs = [
    (0.03, 0., 2.12, 4., 0),
    (0.047, 0., 2.12, 4., 0),
    (0.047, 0.18, 2.12, 4., 0),
    (0.047, 0.18, 1.27, 4., 0),
    (0.047, 0.18, 1.27, 6.9, 0),
]

Qlog0 = [
    (0.03, 0., 2.12, 0.03, 0),
    (0.047, 0., 2.12, 0.047, 0),
    (0.047, 0.18, 2.12, 0.047, 0),
    (0.047, 0.18, 1.27, 0.047, 0),
]

Qlog2 = [
    (0.03, 0.001, 2.12, 4., 2.),
    (0.047, 0.001, 2.12, 4., 2.),
    (0.047, 0.18, 2.12, 4., 2.),
    (0.047, 0.18, 1.27, 4., 2.),
    (0.047, 0.18, 1.27, 6.9, 2.),
    (0.047, 0.18, 1.27, 6.9, 3.21),
]

Qlogs2 = [
    (0.03, 0.001, 2.12, 4., 0),
    (0.047, 0.001, 2.12, 4., 0),
    (0.047, 0.18, 2.12, 4., 0),
    (0.047, 0.18, 1.27, 4., 0),
    (0.047, 0.18, 1.27, 6.9, 0),
]

Qlog02 = [
    (0.03, 0.001, 2.12, 0.03, 0),
    (0.047, 0.001, 2.12, 0.047, 0),
    (0.047, 0.18, 2.12, 0.047, 0),
    (0.047, 0.18, 1.27, 0.047, 0),
]

Q = [(0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
     (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
     (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
     (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
     (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
     (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
     (0.35, 0.75, 1.15, 3.76, 2.22, 0.19)]

Qsm = [(0.62, 0.21, 2.16, 1.2, 1.2, 3.2),
       (0.35, 0.21, 2.16, 1.2, 1.2, 3.2),
       (0.35, 0.75, 2.16, 1.2, 1.2, 3.2),
       (0.35, 0.75, 1.15, 1.2, 1.2, 3.2),
       (0.35, 0.75, 1.15, 3.76, 3.76, 3.2),
       (0.35, 0.75, 1.15, 3.76, 3.76, 0.19)]

Qm0 = [(0.62, 0.21, 2.16, 1.2, 0, 3.2),
       (0.35, 0.21, 2.16, 1.2, 0, 3.2),
       (0.35, 0.75, 2.16, 1.2, 0, 3.2),
       (0.35, 0.75, 1.15, 1.2, 0, 3.2),
       (0.35, 0.75, 1.15, 3.76, 0, 3.2),
       (0.35, 0.75, 1.15, 3.76, 0, 0.19)]

Q00 = [
    (0.62, 0.21, 2.16, 0., 0., 3.2),
    (0.35, 0.21, 2.16, 0., 0., 3.2),
    (0.35, 0.75, 2.16, 0., 0., 3.2),
    (0.35, 0.75, 1.15, 0., 0., 3.2),
    (0.35, 0.75, 1.15, 0., 0., 0.19),
]

QDJm = [
    (0.35, 0.1, 1.2),
    (0.72, 0.1, 1.2),
    (0.72, 0.42, 1.2),
    (0.72, 0.42, 3.18),
]

QDJ0 = [
    (0.35, 0.0, 1.2),
    (0.72, 0.0, 1.2),
    (0.72, 0.0, 3.18),
]

Qlogzm = [
    (0.03, 0, 2.12, 0, 2.),
    (0.047, 0, 2.12, 0, 2.),
    (0.047, 0, 1.27, 0, 2.),
    (0.047, 0, 1.27, 0, 3.21),
]

Qlogzmsm = [(0.03, 0, 2.12, 0, 0),
            (0.047, 0, 2.12, 0, 0),
            (0.047, 0, 1.27, 0, 0)]

Qlogzm0 = [(0.03, 0, 2.12, 0, 0),
           (0.047, 0, 2.12, 0, 0),
           (0.047, 0, 1.27, 0, 0)]

Qzm = [
    (0.62, 0, 2.16, 1.2, 0.8, 3.2),
    (0.35, 0, 2.16, 1.2, 0.8, 3.2),
    (0.35, 0, 1.15, 1.2, 0.8, 3.2),
    (0.35, 0, 1.15, 3.76, 0.8, 3.2),
    (0.35, 0, 1.15, 3.76, 2.22, 3.2),
    (0.35, 0, 1.15, 3.76, 2.22, 0.19),
]

Qzmsm = [
    (0.62, 0, 2.16, 1.2, 1.2, 3.2),
    (0.35, 0, 2.16, 1.2, 1.2, 3.2),
    (0.35, 0, 1.15, 1.2, 1.2, 3.2),
    (0.35, 0, 1.15, 3.76, 3.76, 3.2),
    (0.35, 0, 1.15, 3.76, 3.76, 0.19),
]

Qzmm0 = [
    (0.62, 0, 2.16, 1.2, 0, 3.2),
    (0.35, 0, 2.16, 1.2, 0, 3.2),
    (0.35, 0, 1.15, 1.2, 0, 3.2),
    (0.35, 0, 1.15, 3.76, 0, 3.2),
    (0.35, 0, 1.15, 3.76, 0, 0.19),
]

Qzm00 = [
    (0.62, 0, 2.16, 0, 0, 3.2),
    (0.35, 0, 2.16, 0, 0, 3.2),
    (0.35, 0, 1.15, 0, 0, 3.2),
    (0.35, 0, 1.15, 0, 0, 0.19),
]

Qzp = [(0.62, 0.21, 1.2, 0.8, 3.2),
       (0.35, 0.21, 1.2, 0.8, 3.2),
       (0.35, 0.75, 1.2, 0.8, 3.2),
       (0.35, 0.75, 3.76, 0.8, 3.2),
       (0.35, 0.75, 3.76, 2.22, 3.2),
       (0.35, 0.75, 3.76, 2.22, 0.19)]

Qzpsm = [(0.62, 0.21,  1.2,  1.2, 3.2),
         (0.35, 0.21,  1.2,  1.2, 3.2),
         (0.35, 0.75,  1.2,  1.2, 3.2),
         (0.35, 0.75,  3.76, 3.76, 3.2),
         (0.35, 0.75,  3.76, 3.76, 0.19)]

Qzpm0 = [(0.62, 0.21, 1.2,  0, 3.2),
         (0.35, 0.21, 1.2,  0, 3.2),
         (0.35, 0.75, 1.2,  0, 3.2),
         (0.35, 0.75, 3.76, 0, 3.2),
         (0.35, 0.75, 3.76, 0, 0.19)]

Qzp00 = [(0.62, 0.21, 0, 0, 3.2),
         (0.35, 0.21, 0, 0, 3.2),
         (0.35, 0.75, 0, 0, 3.2),
         (0.35, 0.75, 0, 0, 0.19)]

Qig = [(0.62, 0.21, 2.16, 1.2, 3.2),
       (0.35, 0.21, 2.16, 1.2, 3.2),
       (0.35, 0.75, 2.16, 1.2, 3.2),
       (0.35, 0.75, 1.15, 1.2, 3.2),
       (0.35, 0.75, 1.15, 3.76, 3.2),
       (0.35, 0.75, 1.15, 3.76, 0.19)]

Qigzm = [(0.62, 0, 2.16, 1.2, 3.2),
         (0.35, 0, 2.16, 1.2, 3.2),
         (0.35, 0, 1.15, 1.2, 3.2),
         (0.35, 0, 1.15, 3.76, 3.2),
         (0.35, 0, 1.15, 3.76, 0.19)]

Qigzp = [(0.62, 0.21, 1.2, 3.2),
         (0.35, 0.21, 1.2, 3.2),
         (0.35, 0.75, 1.2, 3.2),
         (0.35, 0.75, 3.76, 3.2),
         (0.35, 0.75, 3.76, 0.19)]

Qigq = [(0.62, 0.21, 2.16, 1.2, 3.2, 0.8),
        (0.35, 0.21, 2.16, 1.2, 3.2, 0.8),
        (0.35, 0.75, 2.16, 1.2, 3.2, 0.8),
        (0.35, 0.75, 1.15, 1.2, 3.2, 0.8),
        (0.35, 0.75, 1.15, 3.76, 3.2, 0.8),
        (0.35, 0.75, 1.15, 3.76, 0.19, 0.8),
        (0.35, 0.75, 1.15, 3.76, 0.19, 5.9)]

Qigqzm = [(0.62, 0, 2.16, 1.2, 3.2, 0.8),
          (0.35, 0, 2.16, 1.2, 3.2, 0.8),
          (0.35, 0, 1.15, 1.2, 3.2, 0.8),
          (0.35, 0, 1.15, 3.76, 3.2, 0.8),
          (0.35, 0, 1.15, 3.76, 0.19, 0.8),
          (0.35, 0, 1.15, 3.76, 0.19, 5.9)]

Qigqzp = [(0.62, 0.21, 1.2, 3.2, 0.8),
          (0.35, 0.21, 1.2, 3.2, 0.8),
          (0.35, 0.75, 1.2, 3.2, 0.8),
          (0.35, 0.75, 3.76, 3.2, 0.8),
          (0.35, 0.75, 3.76, 0.19, 0.8),
          (0.35, 0.75, 3.76, 0.19, 5.9)]

# reslog = [tlog(*args) for args in Qlog]
# reslogs = [tlog(*args) for args in Qlogs]
# reslog0 = [tlog(*args) for args in Qlog0]
# res = [i_m_m_i(*args) for args in Q]
# ressm = [i_m_m_i(*args) for args in Qsm]
# resm0 = [i_m_m_i(*args) for args in Qm0]
# res00 = [i_m_m_i(*args) for args in Q00]
# reslog_l = [tlog_l2(*args) for args in Qlog2]
# reslogs_l = [tlog_l2(*args) for args in Qlogs2]
# reslog0_l = [tlog_l2(*args) for args in Qlog02]
# resl = [i_m_m_l_i(*args) for args in Q]
# reslsm = [i_m_m_l_i(*args) for args in Qsm]
# reslm0 = [i_m_m_l_i(*args) for args in Qm0]
# resl00 = [i_m_m_l_i(*args) for args in Q00]
# reslog_t = [tlog_t(*args) for args in Qlog2]
# reslogs_t = [tlog_t(*args) for args in Qlogs2]
# reslog0_t = [tlog_t(*args) for args in Qlog02]
# rest = [i_m_m_t_i(*args) for args in Q]
# restsm = [i_m_m_t_i(*args) for args in Qsm]
# restm0 = [i_m_m_t_i(*args) for args in Qm0]
# rest00 = [i_m_m_t_i(*args) for args in Q00]

# resdjm = [d_j_m_i(*args) for args in QDJm]
# resdj0 = [d_j_m_i(*args) for args in QDJ0]
# resd2jm = [d2_j_m_i(*args) for args in QDJm]
# resd2j0 = [d2_j_m_i(*args) for args in QDJ0]
# resdjml = [d_j_m_l_i(*args) for args in QDJm]
# resdj0l = [d_j_m_l_i(*args) for args in QDJ0]
# resd2jml = [d2_j_m_l_i(*args) for args in QDJm]
# resd2j0l = [d2_j_m_l_i(*args) for args in QDJ0]
# resdjmt = [d_j_m_t_i(*args) for args in QDJm]
# resdj0t = [d_j_m_t_i(*args) for args in QDJ0]
# resd2jmt = [d2_j_m_t_i(*args) for args in QDJm]
# resd2j0t = [d2_j_m_t_i(*args) for args in QDJ0]

# resd = [d_i_m_m_i(*args) for args in Q]
# resdsm = [d_i_m_m_i(*args) for args in Qsm]
# resdm0 = [d_i_m_m_i(*args) for args in Qm0]
# resd00 = [d_i_m_m_i(*args) for args in Q00]

# resdl = [d_i_m_m_l_i(*args) for args in Q]
# resdlsm = [d_i_m_m_l_i(*args) for args in Qsm]
# resdlm0 = [d_i_m_m_l_i(*args) for args in Qm0]
# resdl00 = [d_i_m_m_l_i(*args) for args in Q00]

# resdt = [d_i_m_m_t_i(*args) for args in Q]
# resdtsm = [d_i_m_m_t_i(*args) for args in Qsm]
# resdtm0 = [d_i_m_m_t_i(*args) for args in Qm0]
# resdt00 = [d_i_m_m_t_i(*args) for args in Q00]

# reslogzm = [tlog(*args) for args in Qlogzm]
# reslogzmsm = [tlog(*args) for args in Qlogzmsm]
# reslogzm0 = [tlog(*args) for args in Qlogzm0]
# reslogtzm = [tlog_t(*args) for args in Qlogzm]
# reslogtzmsm = [tlog_t(*args) for args in Qlogzmsm]
# reslogtzm0 = [tlog_t(*args) for args in Qlogzm0]

# reszm = [i_m_m_i(*args) for args in Qzm]
# reszmsm = [i_m_m_i(*args) for args in Qzmsm]
# reszmm0 = [i_m_m_i(*args) for args in Qzmm0]
# reszm00 = [i_m_m_i(*args) for args in Qzm00]

# reszml = [i_m_m_l_i(*args) for args in Qzm]
# reszmlsm = [i_m_m_l_i(*args) for args in Qzmsm]
# reszmlm0 = [i_m_m_l_i(*args) for args in Qzmm0]
# reszml00 = [i_m_m_l_i(*args) for args in Qzm00]

# reszmt = [i_m_m_t_i(*args) for args in Qzm]
# reszmtsm = [i_m_m_t_i(*args) for args in Qzmsm]
# reszmtm0 = [i_m_m_t_i(*args) for args in Qzmm0]
# reszmt00 = [i_m_m_t_i(*args) for args in Qzm00]

# resdzm = [d_i_m_m_i(*args) for args in Qzm]
# resdzmsm = [d_i_m_m_i(*args) for args in Qzmsm]
# resdzmm0 = [d_i_m_m_i(*args) for args in Qzmm0]
# resdzm00 = [d_i_m_m_i(*args) for args in Qzm00]

# resdlzm = [d_i_m_m_l_i(*args) for args in Qzm]
# resdlzmsm = [d_i_m_m_l_i(*args) for args in Qzmsm]
# resdlzmm0 = [d_i_m_m_l_i(*args) for args in Qzmm0]
# resdlzm00 = [d_i_m_m_l_i(*args) for args in Qzm00]

# resdtzm = [d_i_m_m_t_i(*args) for args in Qzm]
# resdtzmsm = [d_i_m_m_t_i(*args) for args in Qzmsm]
# resdtzmm0 = [d_i_m_m_t_i(*args) for args in Qzmm0]
# resdtzm00 = [d_i_m_m_t_i(*args) for args in Qzm00]

# reszp = [i_m_m_i_zp(*args) for args in Qzp]
# reszpsm = [i_m_m_i_zp(*args) for args in Qzpsm]
# reszpm0 = [i_m_m_i_zp(*args) for args in Qzpm0]
# reszp00 = [i_m_m_i_zp(*args) for args in Qzp00]

# reslzp = [i_m_m_l_i_zp(*args) for args in Qzp]
# reslzpsm = [i_m_m_l_i_zp(*args) for args in Qzpsm]
# reslzpm0 = [i_m_m_l_i_zp(*args) for args in Qzpm0]
# reslzp00 = [i_m_m_l_i_zp(*args) for args in Qzp00]

# restzp = [i_m_m_t_i_zp(*args) for args in Qzp]
# restzpsm = [i_m_m_t_i_zp(*args) for args in Qzpsm]
# restzpm0 = [i_m_m_t_i_zp(*args) for args in Qzpm0]
# restzp00 = [i_m_m_t_i_zp(*args) for args in Qzp00]

# resdzp = [d_i_m_m_i_zp(*args) for args in Qzp]
# resdzpsm = [d_i_m_m_i_zp(*args) for args in Qzpsm]
# resdzpm0 = [d_i_m_m_i_zp(*args) for args in Qzpm0]
# resdzp00 = [d_i_m_m_i_zp(*args) for args in Qzp00]

# resdlzp = [d_i_m_m_l_i_zp(*args) for args in Qzp]
# resdlzpsm = [d_i_m_m_l_i_zp(*args) for args in Qzpsm]
# resdlzpm0 = [d_i_m_m_l_i_zp(*args) for args in Qzpm0]
# resdlzp00 = [d_i_m_m_l_i_zp(*args) for args in Qzp00]

# resghse = [ghost_se(*args) for args in Qig]
# resghsezm = [ghost_se(*args) for args in Qigzm]
# resghsezp = [ghost_se_zp(*args) for args in Qigzp]

# resglpoll = [gluon_pol_l(*args) for args in Qig]
# resglpollzm = [gluon_pol_l(*args) for args in Qigzm]
# resglpollzp = [gluon_pol_l_zp(*args) for args in Qigzp]

# resglpolt = [gluon_pol_t(*args) for args in Qig]
# resglpoltzm = [gluon_pol_t(*args) for args in Qigzm]
# resglpoltzp = [gluon_pol_t_zp(*args) for args in Qigzp]

# resglpolql = [gluon_pol_quark_l(*args) for args in Qigq]
# resglpolqlzm = [gluon_pol_quark_l(*args) for args in Qigqzm]
# resglpolqlzp = [gluon_pol_quark_l_zp(*args) for args in Qigqzp]

# print(str(reslog).replace('j', '*I').replace('(', '').replace(')', ''))
# print(str(reslogs).replace('j', '*I').replace('(', '').replace(')', ''))
# print(str(reslog0).replace('j', '*I').replace('(', '').replace(')', ''))
# print(str(res).replace('j', '*I').replace('(', '').replace(')', ''))
# print(str(ressm).replace('j', '*I').replace('(', '').replace(')', ''))
# print(str(resm0).replace('j', '*I').replace('(', '').replace(')', ''))
# print(str(res00).replace('j', '*I').replace('(', '').replace(')', ''))
# print(str(reslog_l).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslogs_l).replace('j', '*I').replace('(',
#                                                 '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslog0_l).replace('j', '*I').replace('(',
#                                                 '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resl).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslsm).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslm0).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resl00).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslog_t).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslogs_t).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslog0_t).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(rest).replace('j', '*I').replace('(',
#                                            '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(restsm).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(restm0).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(rest00).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdjm).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdj0).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resd2jm).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resd2j0).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdjml).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdj0l).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resd2jml).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resd2j0l).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdjmt).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdj0t).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resd2jmt).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resd2j0t).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resd).replace('j', '*I').replace('(',
#                                           '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdsm).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdm0).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resd00).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdl).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdlsm).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdlm0).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdl00).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdt).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdtsm).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdtm0).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdt00).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslogzm).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslogzmsm).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslogzm0).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslogtzm).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslogtzmsm).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslogtzm0).replace('j', '*I').replace('(',
#      '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszm).replace('j', '*I').replace('(',
#                                            '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszmsm).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszmm0).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszm00).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszml).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszmlsm).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszmlm0).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszml00).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszmt).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszmtsm).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszmtm0).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszmt00).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdzm).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdzmsm).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdzmm0).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdzm00).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdlzm).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdlzmsm).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdlzmm0).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdlzm00).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdtzm).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdtzmsm).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdtzmm0).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdtzm00).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszp).replace('j', '*I').replace('(',
#                                            '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszpsm).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszpm0).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reszp00).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslzp).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslzpsm).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslzpm0).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(reslzp00).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(restzp).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(restzpsm).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(restzpm0).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(restzp00).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I'))
# print(str(resdzp).replace('j', '*I').replace('(',
#                                             '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resdzpsm).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resdzpm0).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resdzp00).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resdlzp).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resdlzpsm).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resdlzpm0).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resdlzp00).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))

# print(str(resghse).replace('j', '*I').replace('(',
#                                              '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resghsezm).replace('j', '*I').replace('(',
#                                               '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resghsezp).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))

# print(str(resglpoll).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resglpollzm).replace('j', '*I').replace('(',
#                                                  '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resglpollzp).replace('j', '*I').replace('(',
#                                                  '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))

# print(str(resglpolt).replace('j', '*I').replace('(',
#                                                '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resglpoltzm).replace('j', '*I').replace('(',
#                                                  '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resglpoltzp).replace('j', '*I').replace('(',
#                                                  '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))

# print(str(resglpolql).replace('j', '*I').replace('(',
#                                                 '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resglpolqlzm).replace('j', '*I').replace('(',
#                                                   '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
# print(str(resglpolqlzp).replace('j', '*I').replace('(',
#                                                   '').replace(')', '').replace('+0*I', '+0.*I').replace('-0*I', '+0.*I'))
