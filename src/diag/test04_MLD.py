import numpy as np
import data_loading_helper as dlh
import MITgcmDiff.loadFunctions as lf
import pandas as pd
import MITgcmDiff.mixed_layer_tools as mlt
from MITgcmDiff.buoyancy_linear import TS2rho
data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
grid_dir = "/data/SO2/SWOT/GRID/BIN"

test_dt    = pd.Timestamp('2017-01-06')
msm = dlh.MITgcmSimMetadata('2015-01-01', 150.0, data_dir, grid_dir)
   
print("Test date: ", test_dt)

lat_rng = [31.0, 43.0]
lon_rng = [230.0, 244.0] #360 .- [130.0, 116.0

#lat_rng = [35.0, 37.0]
##lon_rng = [235.0, 237] #360 .- [130.0, 116.0
nlev = 90

print("Start loading data...")
coo, crop_kwargs = lf.loadCoordinateFromFolderAndWithRange(msm.grid_dir, nlev=nlev, lat_rng=lat_rng, lon_rng=lon_rng)
data = dlh.loadDataByDate(test_dt, msm, **crop_kwargs)
print("Data loaded.")


#print(coo.grid["Depth"])
#print(coo.grid["Nz_bot"])
# Start plotting stuff
#print(coo.grid["RF"])
MLD = mlt.findMLD_drhodz(data["DRHODR"], coo.grid["RF"].flatten(), mask=coo.grid["maskInC"], Nz_bot=coo.grid["Nz_bot"])


rho = TS2rho(data["THETA"], data["SALT"])
MLD_threshold = mlt.findMLD_rho(rho, coo.grid["RC"].flatten(), mask=coo.grid["maskInC"], Nz_bot=coo.grid["Nz_bot"])

lat = coo.grid["YC"][:, 0]
lon = coo.grid["XC"][0, :]

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

fig, ax = plt.subplots(1, 1)

ax.plot(data["DRHODR"][:, 5, 5], coo.grid["RF"].flatten()[:-1])


plt.show()







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
    ncol = 2,
    nrow = 2,
)

fig, ax = plt.subplots(
    2, 2,
    figsize=figsize,
    subplot_kw=dict(projection=proj, aspect="auto"),
    gridspec_kw=gridspec_kw,
    constrained_layout=False,
    squeeze=False,
)

fig.suptitle("Date: %s" % (test_dt.strftime("%Y-%m-%d")))

_ax = ax[0, 0]
mappable = _ax.contourf(lon, lat, coo.grid["Depth"], levels=np.linspace(0, 6000, 601), transform=proj_norm, extend="max", cmap="rainbow")
cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
cb.ax.set_ylabel(" Depth [ $ \\mathrm{m} $ ]")

_ax = ax[0, 1]
mappable = _ax.contourf(lon, lat, coo.grid["Nz_bot"], levels=np.linspace(0, 90, 91), transform=proj_norm, extend="max", cmap="rainbow")
 
cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
cb.ax.set_ylabel(" Nz_bot ")

_ax = ax[1, 0]
mappable = _ax.contourf(lon, lat, MLD, levels=np.linspace(0, 200, 21), transform=proj_norm, extend="max", cmap="rainbow")
 
cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
cb.ax.set_ylabel("MLD [$\\mathrm{m}$]")


_ax = ax[1, 1]
mappable = _ax.contourf(lon, lat, MLD_threshold, levels=np.linspace(0, 200, 21), transform=proj_norm, extend="max", cmap="rainbow")
 
cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
cb.ax.set_ylabel("MLD [$\\mathrm{m}$]")





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
