import numpy as np
import xarray as xr
import pandas as pd

import argparse

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--input', type=str, help='Input file', required=True)
parser.add_argument('--wateryears', type=int, nargs='+', help='Wateryears', required=True)
args = parser.parse_args()
print(args)


ds = xr.open_dataset(args.input)





print("Loading matplotlib...")
import matplotlib.pyplot as plt
from matplotlib import dates
print("Done")


plotting_wateryears = args.wateryears


fig, ax = plt.subplots(3, 1, figsize=(6, 6), gridspec_kw=dict(hspace=.3))

fig.suptitle(args.input)

for wyr in plotting_wateryears:

    print("Plotting: ", wyr)

    time_cond = ( (ds.time.dt.year == wyr-1) & (ds.time.dt.month.isin([10, 11, 12])) ) | ( (ds.time.dt.year == wyr) & (ds.time.dt.month.isin([1, 2, 3])) )


    time_cond = time_cond & np.logical_not( (ds.time.dt.month == 2) & (ds.time.dt.day == 29) )
    
    _ds = ds.where(time_cond, drop=True)

    t = pd.date_range('2001-10-01', periods=len(_ds.time), freq="D")

    ax[0].plot(t, _ds.AR_presense, label="%04d" % (wyr,))
    ax[1].plot(t, _ds.mean_IVT,)
    ax[2].plot(t, _ds.mean_IWV,)

for _ax in ax.flatten():
    _ax.xaxis.set_major_formatter(dates.DateFormatter('%m/%d'))

ax[0].legend()

ax[0].set_title("AR presense (fraction of area that is in the AR objects)")
ax[1].set_title("area mean IVT")
ax[2].set_title("area mean IWV")


fig.savefig("fig/AR_presense.png", dpi=200)


plt.show() 

