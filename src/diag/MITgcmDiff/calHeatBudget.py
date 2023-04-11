import numpy as np
import MITgcmDiff.Operators as op
import MITgcmDiff.xarrayCoversion as xac

rhoConst = 1027.5
c_p      = 3994.0


def Qsw_shape(z):
    return ( 0.62 * np.exp(z/0.6) + (1 - 0.62) * np.exp(z/20.0) ) * (z >= -200.0)

def SFCFLX_shape(z):
    return z >= 0


def cvt2Dto3D(a, coo):
    return np.expand_dims(a, 0)

def computeTendency(
    d,
    coo,
    return_xarray = False,
):
    _d = dict()

    _d["TOTTTEND"] = d["TOTTTEND"] / 86400.0

    _d["SWFLX"]  = - cvt2Dto3D(d["oceQsw"], coo)              * Qsw_shape(coo.grid["RF"][:-1])
    _d["SFCFLX"] = - cvt2Dto3D(d["TFLUX"] - d["oceQsw"], coo) * SFCFLX_shape(coo.grid["RF"][:-1])
    _d["WTHMASS_masked"] = d["WTHMASS"] * SFCFLX_shape(coo.grid["RF"][:-1])

    

    print("Shape of mask")
    print(_d["WTHMASS_masked"].shape)

    #print(Qsw_shape(coo.grid["RF"][:-1])) 
    print("Compute ADVx_TH")
    #print(d["ADVx_TH"].shape)
    #print(coo.grid["DVOLT"].shape)
    _d["TEND_ADVx"] = - (
            op.T_DIVx_U(d["ADVx_TH"], coo, weighted=True)
    )

    print("Compute ADVy_TH")
    _d["TEND_ADVy"] = - (
            op.T_DIVy_V(d["ADVy_TH"], coo, weighted=True)
    )

    print("Compute ADVz_TH")
    _d["TEND_ADVz"] = - (
            op.T_DIVz_W(d["ADVr_TH"], coo, weighted=True)
    )


    _d["TEND_ADV"] = _d["TEND_ADVx"] + _d["TEND_ADVy"] + _d["TEND_ADVz"]


    _d["TEND_DIFFz"] = - (
           op.T_DIVz_W(d["DFrI_TH"], coo, weighted=True)
    )

    _d["TEND_SWFLX"]  = - op.T_DIVz_W(_d["SWFLX"], coo, weighted=False)  / (rhoConst * c_p)
    _d["TEND_SFCFLX"] = - op.T_DIVz_W(_d["SFCFLX"], coo, weighted=False) / (rhoConst * c_p)
    _d["TEND_SFC_WTHMASS"] = - op.T_DIVz_W(_d["WTHMASS_masked"], coo, weighted=False)




    _d["TEND_SUM"] = (
        _d["TEND_ADVx"]
         + _d["TEND_ADVy"]
         + _d["TEND_ADVz"]
         + _d["TEND_DIFFz"]
         + _d["TEND_SWFLX"]
         + _d["TEND_SFCFLX"]
         + _d["TEND_SFC_WTHMASS"]
    )

    _d["TEND_RES"] = _d["TEND_SUM"] - _d["TOTTTEND"]

    if return_xarray:
        
        dict_dimarr = {
                varname : ("T", vardata) for varname, vardata in _d.items()
        }


        ds = xac.mkDataset(
            dict_dimarr = dict_dimarr,
            coo = coo,
        )

        return ds

    else:
        return _d




