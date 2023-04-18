import MITgcmDiff.loadFunctions
import MITgcmDiff.Operators as op 
import xarray as xr
import MITgcmDiff.utils as ut
from MITgcmutils import mds
import MITgcmDiff.loadFunctions as lf
import MITgcmDiff.xarrayCoversion as xac
import MITgcmDiff.calBudget as cb
import numpy as np

import data_loading_helper as dlh


import pandas as pd

data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
grid_dir = "/data/SO2/SWOT/GRID/BIN"

test_dt    = pd.Timestamp('2017-01-06')
msm = dlh.MITgcmSimMetadata('2015-01-01', 150.0, data_dir, grid_dir)
   
print("Test date: ", test_dt)


output_file = "heat_budget_closure.nc"

nlev = 90
lat_rng = [31.0, 43.0]
lon_rng = [230.0, 244.0] #360 .- [130.0, 116.0
#lon_rng = [235.0, 236.0] #360 .- [130.0, 116.0

print("nlev : %d" % (nlev,))
print("Lat range: ", lat_rng)
print("Lon range: ", lon_rng)

print("Loading data")
coo, crop_kwargs = MITgcmDiff.loadFunctions.loadCoordinateFromFolderAndWithRange(grid_dir, nlev=nlev, lat_rng=lat_rng, lon_rng=lon_rng)

data = dlh.loadDataByDate(test_dt, msm, **crop_kwargs, datasets=["diag_Tbdgt", "diag_2D", "diag_state",])

print("Compute heat budget")
ds = cb.computeHeatTendency(data, coo, return_xarray = True)

print("Output file: ", output_file)
ds.to_netcdf(output_file)
