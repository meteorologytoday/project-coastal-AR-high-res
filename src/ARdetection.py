import numpy as np
import xarray as xr
from scipy.ndimage import label, generate_binary_structure
from scipy import spatial


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

            dist_mat[i, j] = getDistOnSphere(candidates[i, 0], candidates[i, 1], candidates[j, 0], candidates[j, 1])

    # get indices of candidates that are furthest apart
    i, j = np.unravel_index(dist_mat.argmax(), dist_mat.shape)
            
    farthest_pair = ( candidates[i, :], candidates[j, :] )

    return farthest_pair, dist_mat[i, j]


def detectARObjects(IVT, coord_lat, coord_lon, area, IVT_threshold):
 
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

        

        centroid = (
            np.sum(coord_lat[idx] * covered_area) / sum_covered_area,
            np.sum(coord_lon[idx] * covered_area) / sum_covered_area,
        )


        AR_objs.append(dict(
            feature_n     = feature_n,
            area          = sum_covered_area,
            centroid      = centroid,
            length        = farthest_dist,
            farthest_pair = farthest_pair,
        ))
    
    
    return labeled_array, AR_objs



# Algorithm


if __name__  == "__main__" :
    
    import xarray as xr

    test_file = "./data/ERA5/AR_processed/ERA5_AR_2018-12-24.nc"
    
    ds = xr.open_dataset(test_file)
   
    coord_lat, coord_lon = np.meshgrid(ds.coords["lat"], ds.coords["lon"], indexing='ij')

    coord_lon = coord_lon % 360.0



    dlat = np.deg2rad((ds.coords["lat"][0] - ds.coords["lat"][1]).to_numpy())
    dlon = np.deg2rad((ds.coords["lon"][1] - ds.coords["lon"][0]).to_numpy())

    R_earth = 6.4e6
 
    area = R_earth**2 * np.cos(np.deg2rad(coord_lat)) * dlon * dlat

    print(area)
   
    print("Compute AR_objets") 
    labeled_array, AR_objs = detectARObjects(ds.IVT[0, :, :], coord_lat, coord_lon, area, IVT_threshold=250.0)

    for i, AR_obj in enumerate(AR_objs):
        pts = AR_obj["farthest_pair"]
        cent = AR_obj["centroid"]
        print("[%i] cent=(%f, %f), area=%fkm^2" % (i, cent[0], cent[1], AR_obj["area"] / 1e6))


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

    cent_lon = 0.0

    plot_lon_l = -180.0
    plot_lon_r = 180.0
    plot_lat_b = 0.0
    plot_lat_t = 70.0

    proj = ccrs.PlateCarree(central_longitude=cent_lon)
    proj_norm = ccrs.PlateCarree()

    fig, ax = plt.subplots(
        1, 1,
        figsize=(6, 4),
        subplot_kw=dict(projection=proj),
        gridspec_kw=dict(hspace=0, wspace=0.2),
        constrained_layout=False,
    )

    cmap = cm.get_cmap("ocean_r")

    labeled_array = labeled_array.astype(float)
    labeled_array[labeled_array!=0] = 1.0
    labeled_array[labeled_array==0] = np.nan


    mappable = ax.contourf(ds.coords["lon"], ds.coords["lat"], labeled_array, cmap=cmap,  transform=proj_norm)

    plt.colorbar(mappable, ax=ax, orientation="vertical")

    for i, AR_obj in enumerate(AR_objs):
        pts = AR_obj["farthest_pair"]
        cent = AR_obj["centroid"]
        ax.plot([pts[0][1], pts[1][1]], [pts[0][0], pts[1][0]], 'r-', transform=ccrs.Geodetic())

        ax.text(cent[1], cent[0], "%d" % (i+1), va="center", ha="center", color="cyan", transform=proj_norm)

    ax.set_global()
    ax.coastlines()
    ax.set_extent([plot_lon_l, plot_lon_r, plot_lat_b, plot_lat_t], crs=proj_norm)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
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


