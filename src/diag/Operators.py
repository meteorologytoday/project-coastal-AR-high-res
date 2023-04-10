import numpy as np
import mitgcm_helper
from MITgcmCoordinate import *

def shift(a, shift, axis, boundary='periodic', fill_value = 0.0, **kwargs):

    if shift == 0:
        return a.copy()
    
    _a = np.roll(a, shift, axis=axis, **kwargs)
    if boundary == 'periodic':

        pass

    elif boundary == 'fill':
        
        indexing = [ slice(None, None, None) for i in range(len(a.shape)) ]

        if shift >= 0:
            indexing[axis] = slice(0, shift)

        elif shift < 0:
            indexing[axis] = slice(a.shape[axis] - ( - shift ), a.shape[axis])

        print(indexing) 
        _a[tuple(indexing)] = fill_value

    else:
        raise Exception("Unknown boundary condition: %s" % (boundary,))
    
    return _a

# Convention follows
#https://mitgcm.readthedocs.io/en/latest/outp_pkgs/outp_pkgs.html#mitgcm-kernel-available-diagnostics-list

def T_DIVx_U(fi, coo: MITgcmCoordinate, weighted = True, boundary='fill'):

    if weighted:    
        fo = ( shift(fi, -1, axis=2, boundary=boundary) - fi ) / coo.grid["DVOLT"]
    else:
        _fi = fi * coo.grid["DYG"]
        fo = ( shift(_fi, -1, axis=2, boundary=boundary) - _fi ) / coo.grid["RAC_slab"]


    return fo


def T_DIVy_V(fi, coo: MITgcmCoordinate, weighted = True, boundary='fill'):

    if weighted:    
        fo = ( shift(fi, -1, axis=1, boundary=boundary) - fi ) / coo.grid["DVOLT"]
    else:
        _fi = fi * coo.grid["DXG"]
        fo = ( shift(_fi, -1, axis=1, boundary=boundary) - _fi ) / coo.grid["RAC_slab"]


    return fo


def T_DIVz_W(fi, coo: MITgcmCoordinate, weighted = True, boundary='fill'):

    if weighted:
        fo = ( shift(fi, -1, axis=0) - fi, boundary=boundary ) / coo.grid["DVOLT"]
    else:
        fo = ( shift(fi, -1, axis=0) - fi, boundary=boundary ) / coo.grid["DRF"]


    return fo


    
 






