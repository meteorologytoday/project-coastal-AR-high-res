
class Coordinate:

    def __init__(self, grid):

        self.grid = grid

        self.Ny = grid["DXG"].shape[0]
        self.Nx = grid["DXG"].shape[1]
        self.Nz = grid["RC"].shape[0]

        #self.grid["DXG_p1"] = np.roll(self.grid["DXG"], -1, axis=0)
        #self.grid["DYG_p1"] = np.roll(self.grid["DYG"], -1, axis=1)

    
        self.grid["RAC_slab"] = self.grid["RAC"].reshape((1, self.Ny, self.Nx))
        self.grid["DVOLT"] = self.grid["DRF"] * self.grid["RAC_slab"]


