from __future__ import division
import const

def cmm1_to_ev(ev):
    return ev*const.hh*const.cc/const.ee*1e2

def ev_to_cmm1(cmm1):
    return cmm1*const.ee/(const.hh*const.cc)*1e-2

def ev_to_wl(ev):
    return const.hh*const.cc/(const.ee*ev)

def wl_to_ev(wl):
    return (const.hh*const.cc)/(const.ee*wl)
