import MITgcmDiff.loadFunctions
import MITgcmDiff.Operators as op 
import xarray as xr
import MITgcmDiff.utils as ut
from MITgcmutils import mds
import MITgcmDiff.loadFunctions as lf
import MITgcmDiff.xarrayCoversion as xac
import MITgcmDiff.calHeatBudget as chb
import numpy as np

output_file = "heat_budget_closure.nc"

data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/TEST_TFLUX/" ; iters=150336;
grid_dir = "/data/SO2/SWOT/GRID/BIN"

nlev = 90
lev = list(range(nlev))

lat_rng = [31.0, 43.0]
lon_rng = [230.0, 244.0] #360 .- [130.0, 116.0
#lon_rng = [235.0, 236.0] #360 .- [130.0, 116.0

print("nlev : %d" % (nlev,))
print("Lat range: ", lat_rng)
print("Lon range: ", lon_rng)

print("Loading coordinate")
coo, crop_kwargs = MITgcmDiff.loadFunctions.loadCoordinateFromFolderAndWithRange(grid_dir, nlev=nlev, lat_rng=lat_rng, lon_rng=lon_rng)

lat = coo.grid["YC"][:, 0]
lon = coo.grid["XC"][0, :]


print("Loading data")

data = dict()
for k in ["diag_state", "diag_Tbdgt",]:
    print("Loading file of ", k)
    bundle = mds.rdmds("%s/%s" % (data_dir, k,), iters, **crop_kwargs, returnmeta=True)
    _data = lf.postprocessRdmds(bundle)

    for varname, vardata in _data.items():
        data[varname] = vardata

for k in ["diag_2D", ]:
    print("Loading file of ", k)
    bundle = mds.rdmds("%s/%s" % (data_dir, k,), iters, region=crop_kwargs["region"], returnmeta=True)
    _data = lf.postprocessRdmds(bundle)

    for varname, vardata in _data.items():
        data[varname] = vardata



ds = chb.computeTendency(data, coo, return_xarray = True)

print("Output file: ", output_file)
ds.to_netcdf(output_file)
