import numpy as np
import xarray as xr
import argparse
import pandas as pd
from pathlib import Path
import tool_fig_config

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--date-rng', type=str, nargs=2, help='Date range.', required=True)
parser.add_argument('--avg-days', type=int, help='How many days before and after to average.', required=True)
parser.add_argument('--input-dir', type=str, help='Input dir', default="output_diag_budgets")
parser.add_argument('--output-dir', type=str, help='Output dir', default="output_figure")
parser.add_argument('--lat', type=float, nargs=2, help='Latitudes in degree', default=[30, 45])
parser.add_argument('--lon', type=float, nargs=2, help='Longitudes in degree', default=[360-130, 360-120])
parser.add_argument('--no-display', action="store_true")
args = parser.parse_args()
print(args)


dts = pd.date_range(args.date_rng[0], args.date_rng[1], freq="D", inclusive="both")

args.lon = np.array(args.lon) % 360.0

lat_n, lat_s = np.amax(args.lat), np.amin(args.lat)
lon_w, lon_e = np.amin(args.lon), np.amax(args.lon)


print("==================================")
print("Date range: ", dts[0], " to ", dts[-1])
print("Latitude  box: %.2f %.2f" % (lat_s, lat_n))
print("Longitude box: %.2f %.2f" % (lon_w, lon_e))
print("==================================")





import matplotlib as mpl
if args.no_display is False:
    mpl.use('TkAgg')
else:
    mpl.use('Agg')
    #mpl.rc('font', size=20)
 
# This program plots the AR event
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.transforms as transforms
from matplotlib.dates import DateFormatter
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

print("done")

print("Create dir: %s" % (args.output_dir,))
Path(args.output_dir).mkdir(parents=True, exist_ok=True)

G_levs  = np.linspace(-1, 1, 51) * 0.5
G_ticks = np.linspace(-1, 1, 11) * 0.5

dMLDdt_levs  = np.linspace(-1, 1, 51) * 10
dMLDdt_ticks  = np.linspace(-1, 1, 21) * 10

MLD_levs  = np.linspace(0,  1, 51) * 100
MLD_ticks  = np.linspace(0, 1, 11) * 100


for dt in dts:

    print("Doing date: ", dt)

    output_filename = "%s/budget_analysis_avg-%d_%s.png" % (args.output_dir, args.avg_days, dt.strftime("%Y-%m-%d"))
    
    _data = []
    
    # filenames

    filenames = [
        "%s/budget_analysis_%s.nc" % (args.input_dir, (dt + (i - args.avg_days) * pd.Timedelta(days=1)).strftime("%Y-%m-%d")) for i in range(2*args.avg_days+1)
    ]


    
    
    print("List of loaded files: ")
    for i, filename in enumerate(filenames):
        print("[%2d] %s" % (i, filename,))

    try:
        _data = xr.open_mfdataset(filenames)

        _data = xr.merge([
            _data,
            (_data["G_sw"] + _data["G_nsw"]).rename("G_sw_nsw"),
        ])

        _data_mean = _data.mean(dim="time")
    except FileNotFoundError as e:
        print("Some files do not exist.")
        print(str(e))
        
        print("Continue to next one...")
        continue

    print("Data loading complete.")
    cent_lon = 180.0

    plot_lon_l = 360 - 130.0
    plot_lon_r = 360 - 120.0
    plot_lat_b = 30.0
    plot_lat_t = 45.0

    proj = ccrs.PlateCarree(central_longitude=cent_lon)
    proj_norm = ccrs.PlateCarree()

    ncol = 4
    nrow = 3

    figsize, gridspec_kw = tool_fig_config.calFigParams(
        w = 3,
        h = 3,
        wspace = 1.5,
        hspace = 0.5,
        w_left = 1.0,
        w_right = 1.5,
        h_bottom = 1.0,
        h_top = 1.0,
        ncol = ncol,
        nrow = nrow,
    )

    fig, ax = plt.subplots(
        nrow, ncol,
        figsize=figsize,
        subplot_kw=dict(projection=proj, aspect="auto"),
        gridspec_kw=gridspec_kw,
        constrained_layout=False,
        squeeze=False,
    )

    coords = _data.coords
    #cmap = cm.get_cmap("bwr")
    #cmap.set_over("green")
    #cmap.set_under("yellow")

    fig.suptitle("%s ; avg days: %d" % ( dt.strftime("%Y-%m-%d"), 1+2*args.avg_days, ))

    ax_i = 0

    _ax = ax.flatten()[ax_i]; ax_i += 1

    mappable = _ax.contourf(coords["lon"], coords["lat"], _data.isel(time=args.avg_days-1).MLT, levels=np.linspace(10, 20, 101), transform=proj_norm, extend="max", cmap="rainbow")
    
    _ax.set_title("MLT")
    """
    _ax.quiver(coords["lon"], coords["lat"], _data.u10.to_numpy(), _data.v10.to_numpy(), scale=200, transform=proj_norm)
 
    cs = _ax.contourf(coords["lon"], coords["lat"], _data['map'], colors='none', levels=[0, 0.5, np.inf], hatches=[None, "."], transform=proj_norm)

    # Remove the contour lines for hatches 
    for _, collection in enumerate(cs.collections):
        collection.set_edgecolor("red")
    """

    cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
    cb = plt.colorbar(mappable, cax=cax, ticks=np.linspace(10, 22, 13), orientation="vertical", pad=0.0)
    cb.ax.set_ylabel(" MLT [ $ {}^\\circ\\mathrm{C} $ ]")

    for i, varname in enumerate(["dMLTdt", "G_adv_fwf", "G_sw", "G_nsw", "G_vdiff", "G_ent", "G_sw_nsw"]): 
        _ax = ax.flatten()[ax_i]; ax_i += 1
        mappable = _ax.contourf(coords["lon"], coords["lat"], _data_mean[varname] * 1e5, levels=G_levs, transform=proj_norm, extend="both", cmap="bwr")
         
        cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
        cb = plt.colorbar(mappable, cax=cax, ticks=G_ticks, orientation="vertical", pad=0.0)
        cb.ax.set_ylabel(" %s [ $ 1 \\times 10^{-5} \\, \\mathrm{K} / \\mathrm{s} $ ]" % (varname,))

        _ax.set_title(varname)


    # dMLDdt        
    _ax = ax.flatten()[ax_i]; ax_i += 1
        
    _ax.set_title("dMLDdt")
    mappable = _ax.contourf(coords["lon"], coords["lat"], _data_mean["dMLDdt"] * 1e5, levels=dMLDdt_levs, transform=proj_norm, extend="both", cmap="BrBG")
     
    cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
    cb = plt.colorbar(mappable, cax=cax, ticks=dMLDdt_ticks, orientation="vertical", pad=0.0)
    cb.ax.set_ylabel(" dMLDdt [ $ 1 \\times 10^{-5} \\, \\mathrm{m} / \\mathrm{s} $ ]")

     # MLD        
    _ax = ax.flatten()[ax_i]; ax_i += 1
        
    _ax.set_title("MLD")
    mappable = _ax.contourf(coords["lon"], coords["lat"], _data_mean["MLD"], levels=MLD_levs, transform=proj_norm, extend="both", cmap="GnBu")
     
    cax = tool_fig_config.addAxesNextToAxes(fig, _ax, "right", thickness=0.03, spacing=0.05)
    cb = plt.colorbar(mappable, cax=cax, ticks=MLD_ticks, orientation="vertical", pad=0.0)
    cb.ax.set_ylabel(" MLD [ $ \\mathrm{m} $ ]")

   
           
    for __ax in ax.flatten()[:ax_i]: 

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


    for i in range(ax_i, len(ax.flatten())):
        ax.flatten()[i].remove()

    print("Output file: ", output_filename)
    
    fig.savefig(output_filename, dpi=200)
    if not args.no_display:
        plt.show()

    plt.close()
    
