import xarray as xr
import numpy as np
import load_data
import pandas as pd


def computeARPresense(IVT, ARobj_map, area, region_idx):

    if np.sum(region_idx) == 0:

        raise Exception("Error: `region_idx` is empty.")
   
    ARidx = (ARobj_map > 0) & region_idx

    sum_area = np.sum(area[region_idx])
    AR_presense = np.sum(area[ARidx]) / sum_area
    mean_IVT  = np.sum(area[ARidx] * IVT[ARidx]) / sum_area

    return dict(
        AR_presense = AR_presense,
        mean_IVT     = mean_IVT,
    )


output_filename = "test.nc"

lat_n = 43.0
lat_s = 31.0
lon_w = -130.0 % 360
lon_e = -120.0 % 360

dts = pd.date_range('2015-01-01', '2019-12-31', freq="D", inclusive="both")

empty_sample = np.zeros((len(dts),))

ds = xr.Dataset(

    data_vars=dict(
        AR_presense = (["time", ], empty_sample.copy()),
        mean_IVT    = (["time", ], empty_sample.copy()),
    ),

    coords=dict(
        time=dts,
    ),

    attrs=dict(description="AR presense data."),

)

area = None

for i, dt in enumerate(dts):
   
    print("Doing date: %s" % (dt.strftime("%Y-%m-%d"),)) 
    _data = dict()

    for  dsname, varname, savename, opts in [
        ("ERA5", "IVT", None, dict()),
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

    result = computeARPresense(_data["IVT"].to_numpy(), _data["map"].to_numpy(), area, region_idx)


    ds.AR_presense[i] = result["AR_presense"] 
    ds.mean_IVT[i]    = result["mean_IVT"]


ds.to_netcdf(
    output_filename,
    unlimited_dims=["time",],
    encoding={'time': {'dtype': 'i4'}},
)


# Load IVT field
# Load AR obj


# Do AR existence
# Do Mean IVT


# Output data

