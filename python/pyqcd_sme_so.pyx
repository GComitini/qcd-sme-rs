from qcd_sme cimport C as FfiC
from qcd_sme cimport ym__gluon__dressing_inv_landau as ffi__ym__gluon__dressing_inv_landau
from qcd_sme cimport ym__gluon__dressing_inv_landau__complex as ffi__ym__gluon__dressing_inv_landau__complex
from qcd_sme cimport ym__gluon__dressing_landau as ffi__ym__gluon__dressing_landau
from qcd_sme cimport ym__gluon__dressing_landau__complex as ffi__ym__gluon__dressing_landau__complex

def ym__gluon__dressing_inv_landau(s, F0):
    return ffi__ym__gluon__dressing_inv_landau(s, F0)

def ym__gluon__dressing_inv_landau__complex(s, F0):
    res = ffi__ym__gluon__dressing_inv_landau__complex(FfiC(s.real, s.imag), F0)
    return res.re+1j*res.im

def ym__gluon__dressing_landau(s, F0):
    return ffi__ym__gluon__dressing_landau(s, F0)

def ym__gluon__dressing_landau__complex(s, F0):
    res = ffi__ym__gluon__dressing_landau__complex(FfiC(s.real, s.imag), F0)
    return res.re+1j*res.im