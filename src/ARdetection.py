import numpy as np
import xarray as xr
from scipy.ndimage import label, generate_binary_structure
from scipy import spatial

r_earth = 6.371e6


def getDistOnSphere(lat1, lon1, lat2, lon2, r=1.0):

    _lat1 = np.deg2rad(lat1)
    _lat2 = np.deg2rad(lat2)
    
    _lon1 = np.deg2rad(lon1)
    _lon2 = np.deg2rad(lon2)

    cosine = (
        np.cos(_lat1) * np.cos(_lat2) * np.cos(_lon1 - _lon2)
        + np.sin(_lat1) * np.sin(_lat2)
    )

    arc = np.arccos(cosine)

    return r * arc
    
    


"""
    `pts` must have the shape of (npts, dim), where dim=2 in AR detection
    
    This algorithm is copied from 
    https://stackoverflow.com/questions/50468643/finding-two-most-far-away-points-in-plot-with-many-points-in-python
"""
def getTheFarthestPtsOnSphere(pts):

    # Looking for the most distant points
    # two points which are fruthest apart will occur as vertices of the convex hulil
    candidates = pts[spatial.ConvexHull(pts).vertices, :]
    #candidates = pts

    # get distances between each pair of candidate points
    # dist_mat = spatial.distance_matrix(candidates, candidates)

    dist_mat = np.zeros((len(candidates), len(candidates)))

    for i in range(len(candidates)):
        for j in range(len(candidates)):

            if i >= j:
                dist_mat[i, j] = 0.0
                continue

            dist_mat[i, j] = getDistOnSphere(candidates[i, 0], candidates[i, 1], candidates[j, 0], candidates[j, 1], r=r_earth)

    # get indices of candidates that are furthest apart
    i, j = np.unravel_index(dist_mat.argmax(), dist_mat.shape)
            
    farthest_pair = ( candidates[i, :], candidates[j, :] )

    return farthest_pair, dist_mat[i, j]


def detectARObjects(IVT, coord_lat, coord_lon, area, IVT_threshold, weight=None, filter_func=None):
 
    # 1. Generate object maps
    # 2. Compute objects' characteristics
       
    IVT_binary = np.zeros(IVT.shape, dtype=int)
    IVT_binary[IVT >= IVT_threshold] = 1    
   
    # Using the default connectedness: four sides
    labeled_array, num_features = label(IVT_binary)

    AR_objs = []


    for feature_n in range(1, num_features+1): # numbering starts at 1 
        
        idx = labeled_array == feature_n
        covered_area = area[idx]
        sum_covered_area = np.sum(covered_area)
        
        Npts = np.sum(idx)
        pts = np.zeros((np.sum(idx), 2))
        pts[:, 0] = coord_lat[idx]
        pts[:, 1] = coord_lon[idx]
        
        if Npts >= 10:
            
            farthest_pair, farthest_dist = getTheFarthestPtsOnSphere(pts)
            
        else:
           
            farthest_dist = 0.0 
            farthest_pair = ( pts[0, :], pts[0, :] )

        if weight is None:
            _wgt = covered_area

        else:
            _wgt = covered_area * weight[idx]

        _sum_wgt = np.sum(_wgt)
        
            
             

        centroid = (
            np.sum(coord_lat[idx] * _wgt) / _sum_wgt,
            np.sum(coord_lon[idx] * _wgt) / _sum_wgt,
        )

        AR_obj = dict(
            feature_n     = feature_n,
            area          = sum_covered_area,
            centroid      = centroid,
            length        = farthest_dist,
            farthest_pair = farthest_pair,
        )
 
        if (filter_func is not None) and (filter_func(AR_obj) is False):
            labeled_array[labeled_array == feature_n] = 0.0
            continue 

        AR_objs.append(AR_obj)
    
    
    return labeled_array, AR_objs


def basicARFilter(AR_obj):

    result = True

    if AR_obj['length'] < 1000e3:
        
        result = False
    
    return result

# Algorithm


if __name__  == "__main__" :
    
    import xarray as xr

    test_file = "./data/ERA5/AR_processed/ERA5_AR_2016-01-15.nc"
    test_clim_file = "./data/ERA5/AR_processed_clim/ERA5_AR_01-15.nc"
    
    ds = xr.open_dataset(test_file)
    ds_clim = xr.open_dataset(test_clim_file)

    print(ds)
    print(ds_clim)

    # find the lon=0
    lon_first_zero = np.argmax(ds.coords["lon"].to_numpy() >= 0)
    print("First longitude zero idx: ", lon_first_zero)
    ds = ds.roll(lon=-lon_first_zero, roll_coords=True)
    ds_clim = ds_clim.roll(lon=-lon_first_zero, roll_coords=True)
    
    lat = ds.coords["lat"].to_numpy() 
    lon = ds.coords["lon"].to_numpy()  % 360
  
    # For some reason we need to reassign it otherwise the contourf will be broken... ??? 
    ds = ds.assign_coords(lon=lon) 
    ds_clim = ds_clim.assign_coords(lon=lon) 
    
    IVT_anom = (ds.IVT - ds_clim.IVT)[0, :, :].to_numpy()
    IVT_full = ds.IVT[0, :, :].to_numpy()

    llat, llon = np.meshgrid(lat, lon, indexing='ij')


    dlat = np.deg2rad((lat[0] - lat[1]))
    dlon = np.deg2rad((lon[1] - lon[0]))

    R_earth = 6.4e6
 
    area = R_earth**2 * np.cos(np.deg2rad(llat)) * dlon * dlat

    print("Compute AR_objets") 
    labeled_array, AR_objs = detectARObjects(IVT_anom, llat, llon, area, IVT_threshold=250.0, weight=IVT_full, filter_func = basicARFilter)
    
    labeled_array_full, AR_objs_full = detectARObjects(IVT_full, llat, llon, area, IVT_threshold=500.0, weight=IVT_full, filter_func = basicARFilter)

    """
    for i, AR_obj in enumerate(AR_objs):
        pts = AR_obj["farthest_pair"]
        cent = AR_obj["centroid"]
        print("[%i] cent=(%f, %f), area=%fkm^2" % (i, cent[0], cent[1], AR_obj["area"] / 1e6))
    """

    print("labeled_array: ", labeled_array.shape) 
    print("Loading matplotlib") 
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.patches import Rectangle
    import matplotlib.transforms as transforms
    from matplotlib.dates import DateFormatter
    import matplotlib.ticker as mticker
    import cartopy.crs as ccrs
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    print("done")

    cent_lon = 180.0

    plot_lon_l = -180.0
    plot_lon_r = 180.0
    plot_lat_b = 0.0
    plot_lat_t = 70.0

    proj = ccrs.PlateCarree(central_longitude=cent_lon)
    proj_norm = ccrs.PlateCarree()

    fig, ax = plt.subplots(
        2, 1,
        figsize=(6, 4),
        subplot_kw=dict(projection=proj),
        gridspec_kw=dict(hspace=0, wspace=0.2),
        constrained_layout=False,
        squeeze=False,
    )



        
    labeled_array = labeled_array.astype(float)
    labeled_array_full = labeled_array_full.astype(float)
        
    labeled_array[labeled_array!=0] = 1.0
    labeled_array_full[labeled_array_full!=0] = 1.0


    for i, (_IVT, _AR_objs, _labeled_array) in enumerate([
        (IVT_full, AR_objs_full, labeled_array_full),
        (IVT_anom, AR_objs, labeled_array),
    ]):

        _ax = ax[i, 0]

        _ax.set_title(["FULL", "ANOM"][i])
        
        levs = [np.linspace(0, 1000, 11), np.linspace(-800, 800, 17)][i]
        cmap = cm.get_cmap([ "ocean_r", "bwr_r" ][i])

        mappable = _ax.contourf(lon, lat, _IVT, levels=levs, cmap=cmap,  transform=proj_norm)
        plt.colorbar(mappable, ax=_ax, orientation="vertical")
        _ax.contour(lon, lat, _labeled_array, levels=[0.5,], colors='yellow',  transform=proj_norm, zorder=98, linewidth=1)


        for i, AR_obj in enumerate(_AR_objs):
            pts = AR_obj["farthest_pair"]
            cent = AR_obj["centroid"]
            _ax.plot([pts[0][1], pts[1][1]], [pts[0][0], pts[1][0]], 'r-', transform=ccrs.Geodetic(), zorder=99)

            _ax.text(cent[1], cent[0], "%d" % (i+1), va="center", ha="center", color="cyan", transform=proj_norm, zorder=100)

        _ax.set_global()
        _ax.coastlines()
        _ax.set_extent([plot_lon_l, plot_lon_r, plot_lat_b, plot_lat_t], crs=proj_norm)

        gl = _ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top   = False
        gl.ylabels_right = False

        #gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
        #gl.xlocator = mticker.FixedLocator([120, 150, 180, -150, -120])#np.arange(-180, 181, 30))
        gl.ylocator = mticker.FixedLocator([10, 20, 30, 40, 50])

        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'black'}
        gl.ylabel_style = {'size': 10, 'color': 'black'}

    plt.show()


