import MITgcmDiff.loadFunctions
import MITgcmDiff.Operators as op 
import xarray as xr
import MITgcmDiff.utils as ut
from MITgcmutils import mds
import MITgcmDiff.loadFunctions as lf
import numpy as np


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


bundle = mds.rdmds("%s/diag_state" % (data_dir,), iters, **crop_kwargs, returnmeta=True)
data = lf.postprocessRdmds(bundle)

print("Loading Matplotlib...")

import matplotlib as mpl
mpl.use('TkAgg')
 
# This program plots the AR event
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.transforms as transforms
from matplotlib.dates import DateFormatter
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import tool_fig_config        

print("done")

cent_lon = 180.0

plot_lon_l = 360 - 130.0
plot_lon_r = 360 - 120.0
plot_lat_b = 30.0
plot_lat_t = 45.0

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()

figsize, gridspec_kw = tool_fig_config.calFigParams(
    w = 5,
    h = 5,
    wspace = 1.5,
    hspace = 0.5,
    w_left = 1.0,
    w_right = 1.5,
    h_bottom = 1.0,
    h_top = 1.0,
    ncol = 1,
    nrow = 1,
)

fig, ax = plt.subplots(
    1, 1,
    figsize=figsize,
    subplot_kw=dict(projection=proj, aspect="auto"),
    gridspec_kw=gridspec_kw,
    constrained_layout=False,
    squeeze=False,
)

_ax = ax[0, 0]
mappable = _ax.contourf(lon, lat, data["THETA"][0, :, :], levels=np.linspace(15, 25, 51), transform=proj_norm, extend="max", cmap="rainbow")
 
cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
cb.ax.set_ylabel(" SST [ $ \\mathrm{C} $ ]")

       
for __ax in ax.flatten(): 

    __ax.set_global()
    #__ax.gridlines()
    __ax.coastlines(color='gray')
    __ax.set_extent([plot_lon_l, plot_lon_r, plot_lat_b, plot_lat_t], crs=proj_norm)

    gl = __ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')

    gl.xlabels_top   = False
    gl.ylabels_right = False

    #gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
    #gl.xlocator = mticker.FixedLocator([120, 150, 180, -150, -120])#np.arange(-180, 181, 30))
    #gl.ylocator = mticker.FixedLocator([10, 20, 30, 40, 50])
    
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

plt.show()
