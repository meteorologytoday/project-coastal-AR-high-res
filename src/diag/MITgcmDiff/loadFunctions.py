from MITgcmutils import mds
import MITgcmDiff.Operators as Operators
from MITgcmDiff.Coordinate import Coordinate

def genDataDictionary(data, metadata):

    
    #    data      :: AbstractArray,
    #    metadata  :: Dict,
    #    copy_data :: Bool = false,
    #)

    dimlist = metadata["dimlist"]
    
    data_dict = dict()
    for (i, fldname) in enumerate(metadata["fldlist"]):

        trimmed_fldname = fldname.lstrip().rstrip()
        data_dict[trimmed_fldname] = data[i]

    return data_dict

 
def postprocessRdmds(bundle):

    data = bundle[0]
    itrs = bundle[1]
    metadata = bundle[2]
    
   
    return genDataDictionary(data, metadata) 


def loadCoordinateFromFolder(dirname, nlev=None, region=None):
    
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

        data[varname] = mds.rdmds("%s/%s" % (dirname, varname,), region=region)


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


    return Coordinate(data)

def loadSkeletonFromFolder(dirname, nlev=None):
    
    data = dict()
    for varname in [
        "XG",
        "YG",
        "XC",
        "YC",
        "RF",
        "RC",
        "maskInC",
        "maskInS",
        "maskInW",
    ]:

        data[varname] = mds.rdmds("%s/%s" % (dirname, varname,))


    if nlev is not None:

        for varname in [
            "RC",
        ]:

            data[varname] = data[varname][0:nlev, :, :]

        for varname in [
            "RF",
        ]:

            data[varname] = data[varname][0:nlev+1, :, :]


    return data




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

    #print(coo.grid["XC"][0, :])
    #print(coo.grid["YC"][:, 0])


    
