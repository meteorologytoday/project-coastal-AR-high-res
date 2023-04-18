import MITgcmDiff.loadFunctions
import MITgcmDiff.Operators as op 
import MITgcmDiff.utils as ut
from MITgcmutils import mds
import MITgcmDiff.loadFunctions as lf
import MITgcmDiff.xarrayCoversion as xac
import numpy as np
import pandas as pd


class MITgcmSimMetadata:

    def __init__(self, start_datetime, deltaT, data_dir, grid_dir):

        self.start_datetime = pd.Timestamp(start_datetime)
        self.deltaT = deltaT
        self.data_dir = data_dir
        self.grid_dir = grid_dir
        

def getItersFromDate(dt, msm : MITgcmSimMetadata):
   
    if type(dt) == str:
        dt = pd.Timestamp(dt) 
 
    delta_seconds = (dt - msm.start_datetime).total_seconds()
    iters = delta_seconds / msm.deltaT

    if iters % 1 != 0:
        raise Exception("The specified time is not a multiple of deltaT. Please check if deltaT is correct.")
    
    iters = int(iters)
    
    return iters

def getDateFromIters(iters, msm : MITgcmSimMetadata):
    
    return msm.start_datetime + pd.Timedelta(seconds=msm.deltaT) * iters




def loadDataByDate(dt, msm : MITgcmSimMetadata, region=None, lev=(), merge=True, datasets=[]):
    
    data = dict()

    iters = getItersFromDate(dt, msm)
    
    for k in datasets:
        
        print("Loading file of ", k)


        kwargs = dict(
            region=region,
            returnmeta=True,
        )

        if k in ["diag_state", "diag_Tbdgt"]:
            kwargs["lev"] = lev

        bundle = mds.rdmds("%s/%s" % (msm.data_dir, k,), iters, **kwargs)
        _data = lf.postprocessRdmds(bundle)


        if merge:
            for varname, vardata in _data.items():
                data[varname] = vardata

        else:
            data[k] = _data


    return data




