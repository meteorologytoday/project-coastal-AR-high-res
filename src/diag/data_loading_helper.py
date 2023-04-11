import MITgcmDiff.loadFunctions
import MITgcmDiff.Operators as op 
import MITgcmDiff.utils as ut
from MITgcmutils import mds
import MITgcmDiff.loadFunctions as lf
import MITgcmDiff.xarrayCoversion as xac
import MITgcmDiff.calHeatBudget as chb
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




def loadSkeleton(msm : MITgcmSimMetadata, nlev : int, lat_rng, lon_rng):
    
    skeleton = MITgcmDiff.loadFunctions.loadSkeletonFromFolder(msm.grid_dir, nlev=nlev)
    
    region = ut.findRegion_latlon(
        skeleton["YC"][:, 0], lat_rng,
        skeleton["XC"][0, :], lon_rng,
    )

    coo = MITgcmDiff.loadFunctions.loadCoordinateFromFolder(msm.grid_dir, nlev=nlev, region=region)

    lev = list(range(nlev))
    
    return coo, dict(region=region, lev=lev)

   
def loadDataByDate(dt, msm : MITgcmSimMetadata, region=None, lev=()):
    
    data = dict()

    iters = getItersFromDate(dt, msm)
    
    for k in ["diag_state", "diag_Tbdgt",]:
        print("Loading file of ", k)
        bundle = mds.rdmds("%s/%s" % (msm.data_dir, k,), iters, region=region, lev=lev, returnmeta=True)
        _data = lf.postprocessRdmds(bundle)

        for varname, vardata in _data.items():
            data[varname] = vardata

    for k in ["diag_2D", ]:
        print("Loading file of ", k)
        bundle = mds.rdmds("%s/%s" % (msm.data_dir, k,), iters, region=region, returnmeta=True)
        _data = lf.postprocessRdmds(bundle)

        for varname, vardata in _data.items():
            data[varname] = vardata

    return data




