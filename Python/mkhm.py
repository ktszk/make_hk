#!/usr/bin/env python
#-*- coding:utf-8 -*-
def mkhm(hop,rvec,klist,nlist):
    import numpy as np
    """
    nr is number of hopping matrix
    no is number of orbitals
    hop: hopping parameters as hop[nr,no,no]
    rvec: space coordinates as rvec[nr,3]
    klist: list of k-point coordinate as [kx,ky,kz]
    nlist: list of number of k-mesh for each axis as [nqx,nqy,nqz]
    """
    pi2=2j*np.pi
    phase=[pi2*sum(r*k/float(n) for r,k,n in zip(rr,klist,nlist)) for rr in rvec] 
    en=np.array([[sum(hp*np.exp(ph)) for hp,ph in zip(hopp,phase)] for hopp in hop])
    return en
