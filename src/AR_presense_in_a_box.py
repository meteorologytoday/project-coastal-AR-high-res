import xarray as xr
import numpy as np
import load_data
import pandas as pd

import argparse

parser = argparse.ArgumentParser(
                    prog = 'plot_skill',
                    description = 'Plot prediction skill of GFS on AR.',
)

parser.add_argument('--lat', type=float, nargs=2, help='Latitudes in degree', required=True)
parser.add_argument('--lon', type=float, nargs=2, help='Longitudes in degree', required=True)
parser.add_argument('--year-rng', type=int, nargs=2, help='Year range.', required=True)
parser.add_argument('--output-file', type=str, help='Output dir', default="")
args = parser.parse_args()
print(args)


def computeARPresense(IVT, IWV, ARobj_map, area, region_idx):

    if np.sum(region_idx) == 0:

        raise Exception("Error: `region_idx` is empty.")
   
    ARidx = (ARobj_map > 0) & region_idx

    sum_area = np.sum(area[region_idx])
    sum_ARobj_area = np.sum(area[ARidx])

    AR_presense = np.sum(area[ARidx]) / sum_area

    mean_IVT = np.sum(area[region_idx] * IVT[region_idx]) / sum_area
    mean_IWV = np.sum(area[region_idx] * IWV[region_idx]) / sum_area

    if sum_ARobj_area == 0:
        mean_IVT_in_ARobj = 0.0
        mean_IWV_in_ARobj = 0.0
    else:
        mean_IVT_in_ARobj = np.sum(area[ARidx] * IVT[ARidx]) / sum_ARobj_area
        mean_IWV_in_ARobj = np.sum(area[ARidx] * IWV[ARidx]) / sum_ARobj_area
    

    return dict(
        AR_presense = AR_presense,
        mean_IVT     = mean_IVT,
        mean_IVT_in_ARobj = mean_IVT_in_ARobj,
        mean_IWV     = mean_IWV,
        mean_IWV_in_ARobj = mean_IWV_in_ARobj,
    )


output_filename = "output.nc"

args.lon = np.array(args.lon) % 360.0


lat_n, lat_s = np.amax(args.lat), np.amin(args.lat)
lon_w, lon_e = np.amin(args.lon), np.amax(args.lon)

print("Latitude  box: %.2f %.2f" % (lat_s, lat_n))
print("Longitude box: %.2f %.2f" % (lon_w, lon_e))


dts = pd.date_range('%04d-10-01' % (args.year_rng[0]-1,), '%04d-03-31' % (args.year_rng[1],), freq="D", inclusive="both")

empty_sample = np.zeros((len(dts),))
empty_sample[:] = np.nan

ds = xr.Dataset(

    data_vars=dict(
        AR_presense = (["time", ], empty_sample.copy()),
        mean_IVT    = (["time", ], empty_sample.copy()),
        mean_IVT_in_ARobj = (["time", ], empty_sample.copy()),
        mean_IWV    = (["time", ], empty_sample.copy()),
        mean_IWV_in_ARobj = (["time", ], empty_sample.copy()),
    ),

    coords=dict(
        time=dts,
    ),

    attrs=dict(description="AR presense data."),

)

area = None

for i, dt in enumerate(dts):
   
    print("Doing date: %s" % (dt.strftime("%Y-%m-%d"),)) 

    if dt.month in [4, 5, 6, 7, 8, 9]:
        print("Skip this date.")
        continue

    _data = dict()

    for  dsname, varname, savename, opts in [
        ("ERA5", "IVT", None, dict()),
        ("ERA5", "IWV", None, dict()),
        ("ERA5_ARobj", "map", None, dict(method="ANOM_LEN")),
    ]:

        
        
        info = load_data.getFileAndIndex(product=dsname, date=dt, varname=varname, **opts)

        if savename is None:
            savename = varname

        _data[savename] = xr.open_dataset(info["filename"])[varname].isel(time=0)


    if area is None:

        __data = _data["IVT"]        
        lat = __data.coords["lat"].to_numpy()
        lon = __data.coords["lon"].to_numpy() % 360.0
        
        llat, llon = np.meshgrid(lat, lon, indexing="ij") 

        dlat = lat[0] - lat[1]  # reversed lat in ERA5 data
        dlon = lon[1] - lon[0]  # reversed lat in ERA5 data

        print("dlat = %.2f, dlon = %.2f" % (dlat, dlon, ))

        area = (6.371e6 ** 2.0) * np.cos(np.deg2rad(llat)) * np.deg2rad(dlat) * np.deg2rad(dlon)



        region_idx = (llat >= lat_s) & (llat <= lat_n) & (llon >= lon_w) & (llon <= lon_e)

    result = computeARPresense(_data["IVT"].to_numpy(), _data["IWV"].to_numpy(), _data["map"].to_numpy(), area, region_idx)


    ds.AR_presense[i]       = result["AR_presense"] 
    ds.mean_IVT[i]          = result["mean_IVT"]
    ds.mean_IVT_in_ARobj[i] = result["mean_IVT_in_ARobj"]
    ds.mean_IWV[i]          = result["mean_IWV"]
    ds.mean_IWV_in_ARobj[i] = result["mean_IWV_in_ARobj"]



print("Output filename: ", args.output_file)


ds.to_netcdf(
    args.output_file,
    unlimited_dims=["time",],
    encoding={'time': {'dtype': 'i4'}},
)

