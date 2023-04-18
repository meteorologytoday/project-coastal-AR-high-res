import numpy as np
import data_loading_helper as dlh
import MITgcmDiff.loadFunctions as lf
import pandas as pd
import MITgcmDiff.mixed_layer_tools as mlt
import MITgcmDiff.calBudget as cb
from MITgcmDiff.buoyancy_linear import TS2rho
data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
grid_dir = "/data/SO2/SWOT/GRID/BIN"

test_dt    = pd.Timestamp('2016-10-14')
msm = dlh.MITgcmSimMetadata('2015-01-01', 150.0, data_dir, grid_dir)
   
print("Test date: ", test_dt)

lat_rng = [31.0, 43.0]
lon_rng = [230.0, 244.0] #360 .- [130.0, 116.0

#lat_rng = [35.0, 37.0]
##lon_rng = [235.0, 237] #360 .- [130.0, 116.0
nlev = 90

print("Start loading data...")
coo, crop_kwargs = lf.loadCoordinateFromFolderAndWithRange(msm.grid_dir, nlev=nlev, lat_rng=lat_rng, lon_rng=lon_rng)

# Load snapshots

data_snp1 = dlh.loadDataByDate(test_dt - pd.Timedelta(days=1), msm, **crop_kwargs, datasets=["diag_snaps",])
data_snp2 = dlh.loadDataByDate(test_dt, msm, **crop_kwargs, datasets=["diag_snaps",])
data_ave  = dlh.loadDataByDate(test_dt, msm, **crop_kwargs, datasets=["diag_Tbdgt", "diag_2D", "diag_state",])

print("Data loaded.")


z_T = coo.grid["RC"].flatten()
z_W = coo.grid["RF"].flatten()
mask = coo.grid["maskInC"]

print("Compute MLD... ")
rho1 = TS2rho(data_snp1["THETA"], data_snp1["SALT"])
rho2 = TS2rho(data_snp2["THETA"], data_snp2["SALT"])
MLD1_threshold = mlt.findMLD_rho(rho1, z_T, mask=mask, Nz_bot=coo.grid["Nz_bot"])
MLD2_threshold = mlt.findMLD_rho(rho2, z_T, mask=mask, Nz_bot=coo.grid["Nz_bot"])


print("Compute heat budget...")
ds = cb.computeHeatTendency(data_ave, coo, return_xarray = True)


print("Compute ML heat budget...")
MLG = dict()
for varnameA, varnameB in [
    ("G_ttl",   "TOTTTEND"),
    ("G_adv",   "TEND_ADV"),
    ("G_vdiff", "TEND_DIFFz"),
    ("G_sw",    "TEND_SWFLX"),
    ("G_nsw",   "TEND_SFCFLX"),
    ("G_fwf",   "TEND_SFC_WTHMASS"),
]:
    print("Var: ", varnameA)
    MLG[varnameA] = mlt.computeMLMean(ds[varnameB].to_numpy(), MLD2_threshold, z_W, mask=mask, fill_value=np.nan)

MLG["dMLTdt"] = (
      mlt.computeMLMean(data_snp2["THETA"], MLD2_threshold, z_W, mask=mask, fill_value=np.nan)  
    - mlt.computeMLMean(data_snp1["THETA"], MLD1_threshold, z_W, mask=mask, fill_value=np.nan)  
) / 86400.0



MLG["G_ent2"] = MLG["dMLTdt"] - ( MLG["G_adv"] + MLG["G_vdiff"] + MLG["G_sw"] + MLG["G_nsw"] + MLG["G_fwf"] )
MLG["G_ent"] = (
      mlt.computeMLMean(data_snp1["THETA"], MLD2_threshold, z_W, mask=mask, fill_value=np.nan)  
    - mlt.computeMLMean(data_snp1["THETA"], MLD1_threshold, z_W, mask=mask, fill_value=np.nan)  
) / 86400.0


MLG["G_ttl_surf"] = ds["TOTTTEND"].to_numpy()[1:2, :, :].mean(axis=0)

dMLDdt = (MLD2_threshold - MLD1_threshold) / 86400.0

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

lat = coo.grid["YC"][:, 0]
lon = coo.grid["XC"][0, :]


cent_lon = 180.0

plot_lon_l = 360 - 130.0
plot_lon_r = 360 - 120.0
plot_lat_b = 30.0
plot_lat_t = 45.0

plot_varnames = ["dMLTdt", "G_adv", "G_vdiff", "G_ent", "G_sw", "G_nsw", "G_fwf", "G_ent2"]

G_lev = np.linspace(-1, 1, 21) * 1e-5

proj = ccrs.PlateCarree(central_longitude=cent_lon)
proj_norm = ccrs.PlateCarree()

figsize, gridspec_kw = tool_fig_config.calFigParams(
    w = 3,
    h = 3,
    wspace = 1.5,
    hspace = 0.5,
    w_left = 1.0,
    w_right = 1.5,
    h_bottom = 1.0,
    h_top = 1.0,
    ncol = 5,
    nrow = 2,
)

fig, ax = plt.subplots(
    2, 5,
    figsize=figsize,
    subplot_kw=dict(projection=proj, aspect="auto"),
    gridspec_kw=gridspec_kw,
    constrained_layout=False,
    squeeze=False,
)

fig.suptitle("Date: %s" % (test_dt.strftime("%Y-%m-%d")))


for i, varname in enumerate(plot_varnames):

    _ax = ax.flatten()[i]
    mappable = _ax.contourf(lon, lat, MLG[varname], levels=G_lev, transform=proj_norm, extend="both", cmap="bwr")
    cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
    cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
    cb.ax.set_ylabel(" %s [ $ \\mathrm{K} / \\mathrm{s} $ ]" % (varname,))


_ax = ax.flatten()[-2]
mappable = _ax.contourf(lon, lat, MLD2_threshold, levels=np.linspace(0, 100, 21), transform=proj_norm, extend="both", cmap="GnBu")
cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
cb.ax.set_ylabel(" MLD [ $ \\mathrm{m} $ ]")

_ax = ax.flatten()[-1]
mappable = _ax.contourf(lon, lat, dMLDdt * 86400, levels=np.linspace(-1, 1, 21) * 50, transform=proj_norm, extend="both", cmap="bwr_r")
cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
cb = plt.colorbar(mappable, cax=cax, orientation="vertical", pad=0.0)
cb.ax.set_ylabel(" dMLDdt [ $ \\mathrm{m} / \\mathrm{day} $ ]")


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
