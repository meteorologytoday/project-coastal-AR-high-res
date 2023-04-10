from MITgcmutils import mds
import Operators
from MITgcmCoordinate import *



def loadCoordinateFromFolder(dirname, nlev=None):
    
    data = dict()
    for varname in [
        "DXG",
        "DYG",
        "DXC",
        "DYC",
        "DRC",
        "DRF",
        "RAC",
        "RAZ",
        "XG",
        "YG",
        "XC",
        "YC",
        "RF",
        "RC",
        "maskInC",
        "maskInS",
        "maskInW",
        "hFacC",
        "hFacW",
    ]:

        data[varname] = mds.rdmds("%s/%s" % (dirname, varname,))


    if nlev is not None:

        for varname in [
            "hFacC",
            "hFacW",
            "RC",
            "DRC",
        ]:

            data[varname] = data[varname][0:nlev, :, :]

        for varname in [
            "RF",
            "DRF",
        ]:

            data[varname] = data[varname][0:nlev+1, :, :]


    return MITgcmCoordinate(data)




if __name__ == "__main__":


    input_dirname = "/data/SO2/SWOT/GRID/BIN"


    coo = loadCoordinateFromFolder(input_dirname, nlev=10)

    for k, v in coo.grid.items():    
        print("%s => " % (k,), v.shape) 
    
    print("Testing shift...")
    for boundary in ["periodic", "fill"]:
        for shift in [0, 1, 5, -1, -5]:
            print("Shift on z axes with shift = %d and boundary = %s" % (shift, boundary, ))
            print(Operators.shift(coo.grid["RF"], shift=shift, axis=0, boundary=boundary))




