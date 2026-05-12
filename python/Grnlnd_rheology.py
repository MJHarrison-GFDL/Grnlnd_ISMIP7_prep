import numpy as np  # numerical library
import xarray as xr  # netCDF library
from pyproj import CRS
from pyproj.enums import WktVersion
from pyproj import Transformer

def find_in_domain(bnds,x,y):
    """
    Find bounding region for x,y grid
    """
    i0=np.where(x>=bnds[0][0])[0][0]
    j0=np.where(y<bnds[0][1])[0][0]
    i1=np.where(x>=bnds[1][0])[0][0]
    j1=np.where(y<bnds[1][1])[0][0]

    return [(i0,j0),(i1,j1)]


def linfit2d(x,y,z):
          """
      Fit a set of points in 2-d to a cartesian plane.

          z = a[0] + a[1]x + a[2]y
          """
          a_init = [0.0,0.0,0.0]
          fitfunc = lambda a, x, y: a[0] + a[1]*x + a[2]*y
          ff= fitfunc(a_init,x,y)
          errfunc = lambda a, x, y, z: z - fitfunc(a,x,y)
          err_res= errfunc(a_init,x,y,z)
          out = optimize.leastsq(errfunc,a_init,args=(x.flatten(),y.flatten(),z.flatten()),full_output=1)
          a_final=out[0]

          return a_final

def coarsen_2d(arr, block_size):
    """
    Coarsens a 2D array by averaging blocks of size block_size.
    Example: block_size=(2, 2) averages 2x2 cells into 1 cell.
    """
    m, n = arr.shape
    by, bx = block_size

    # Reshape into (new_rows, block_height, new_cols, block_width)
    reshaped = arr.reshape(m // by, by, n // bx, bx)

    # Average across the block dimensions (axis 1 and 3)
    return reshaped.mean(axis=(1, 3))







class quadmesh:
    """
    Rectangular grid class
    """
    def __init__(self,name=None,parent=None):
        self.name = name
        if parent==None:
            self.x = sgx1km[:-2] # The parent grid contains an odd number of grid cells which is inconvenient for
                                 # coarsening the data. We are re-defining the grid to instead contain an even number of cells
            self.y = sgy1km[:-2]
            self.level = 1
        else:
            self.x = parent.x[::2]
            self.y = parent.y[::2]
            self.level = parent.level + 1
        self.ni=int(0.5*(self.x.shape[0]-1));self.nj=int(0.5*(self.y.shape[0]-1))

    def print(self):
        print(self.name,' grid dimensions: ',(self.nj,self.ni))
        print(self.x.shape,self.y.shape)
        try:
            print('friction_coefficient : (max/min)',self.fralpha.shape,self.fralpha.max(),self.fralpha.min())
        except:
            pass



    def add_rheology(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.fralpha = fralpha1km[:-1,:-1] # cropping the odd cells off the upper-right side of the array
                print('Added friction coefficient ( for 1km grid mean/std: ',self.fralpha.mean(),self.fralpha.std())
                self.Bbar = Bbar1km[:-1,:-1] # cropping the odd cells off the upper-right side of the array
                print('Added Bbar ( for 1km grid mean/std: ',self.Bbar.mean(),self.Bbar.std())
        else:
            self.fralpha = 0.25*((parent.fralpha[:-1:2,:-1:2]+parent.fralpha[:-1:2,1::2])+\
                             parent.fralpha[1::2,:-1:2]+parent.fralpha[1::2,1::2])
            print('Added fralpha elevation for child grid mean/std: ',self.fralpha.mean(),self.fralpha.std())
            self.Bbar = 0.25*((parent.Bbar[:-1:2,:-1:2]+parent.Bbar[:-1:2,1::2])+\
                             parent.Bbar[1::2,:-1:2]+parent.Bbar[1::2,1::2])
            print('Added Bbar elevation for child grid mean/std: ',self.Bbar.mean(),self.Bbar.std())



    def add_vel(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.vx = vx1km[:-1,:-1] # cropping the odd cells off the upper-right side of the array
                self.vy = vy1km[:-1,:-1] # cropping the odd cells off the upper-right side of the array
                print('Added vx for 1km grid mean/std: ',self.vx.mean(),self.vx.std())
                print('Added vy for 1km grid mean/std: ',self.vy.mean(),self.vy.std())
        else:
            vh = parent.vx*parent.thickness
            self.vx =0.25*((vh[:-1:2,:-1:2]+vh[:-1:2,1::2])+\
                      vh[1::2,:-1:2]+vh[1::2,1::2])/np.maximum(self.thickness,1.e-12)
            vh = parent.vy*parent.thickness
            self.vy =0.25*((vh[:-1:2,:-1:2]+vh[:-1:2,1::2])+\
                      vh[1::2,:-1:2]+vh[1::2,1::2])/np.maximum(self.thickness,1.e-12)
            print('Added vx for child grid mean/std: ',self.vx.mean(),self.vx.std())
            print('Added vy for child grid mean/std: ',self.vy.mean(),self.vy.std())

    def add_mask(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.mask = mask1km[:-1,:-1] # cropping the odd cells off the upper-right side of the array
                print('Added Mask for 1km grid max/min: ',self.mask.max(),self.mask.min())
        else:
            self.mask = np.round(0.25*((parent.mask[:-1:2,:-1:2]+parent.mask[:-1:2,1::2])+\
                             parent.mask[1::2,:-1:2]+parent.mask[1::2,1::2]))
            print('Added mask elevation for child grid max/min: ',self.mask.max(),self.mask.min())

    def to_netcdf(self,path=None):
        xh=self.x[1::2]
        yh=self.y[1::2]
        path_out = self.name+'_rheology.nc'
        if path is not None: path_out=path
        print('saving to netcdf, ',path_out, self.name)
        print(xh.shape)
        print(yh.shape)
        print(self.fralpha.shape)
        dA_fr=xr.DataArray(name='friction_coefficient',data=self.fralpha,dims=["y","x"],coords=dict(y=(["y"],yh),x=(["x"],xh)))
        dA_bb=xr.DataArray(name='Bbar',data=self.Bbar,dims=["y","x"],coords=dict(y=(["y"],yh),x=(["x"],xh)))
        ds_out=xr.merge((dA_fr,dA_bb))
        ds_out.to_netcdf(path_out)

path='INPUT_NS/GreenlandISMIP6_BMa5_Control_drag_weertman_abslog_BMa5_interp_d1.nc'
print('Opening dataset: ',path)
ds=xr.open_dataset(path)
#Cell centers xh,yh (m)
xh=ds['x'].load().data
yh=ds['y'].load().data
ni=xh.shape[0];nj=yh.shape[0]
print('Grid size= ',(ni,nj))
print('x center coord range: ',(xh[0],xh[-1]))
print('y center coord range: ',(yh[0],yh[-1]))
#proj4text=ds['mapping'].proj4text
#proj = CRS.from_proj4(proj4text)
#print(proj.to_wkt(WktVersion.WKT1_GDAL, pretty=True))
# Distance between cell centers dx,dy (m)
dx=xh[1:]-xh[:-1]
dy=yh[1:]-yh[:-1]
# quad cell nodal (q) points
xq1km=xh.copy();yq1km=yh.copy()
dx_2=0.5*dx;dy_2=0.5*dy
xq1km[:-1]=xq1km[:-1]-dx_2;xq1km[-1]=xq1km[-1]-dx_2[-1]
yq1km[:-1]=yq1km[:-1]-dy_2;yq1km[-1]=yq1km[-1]-dy_2[-1]
xq1km=np.concatenate((xq1km,[xq1km[-1]+dx[-1]]))
yq1km=np.concatenate((yq1km,[yq1km[-1]+dy[-1]]))
nip=xq1km.shape[0];njp=yq1km.shape[0]
xbnds = (xq1km[0],xq1km[-1])
ybnds = (yq1km[0],yq1km[-1])
print('Grid cell y coordinate range: ',xbnds)
print('Grid cell y coordinate range: ',ybnds)
# Supergrid (h + q)
sgx1km=np.zeros(2*ni+1);sgx1km[::2]=xq1km;sgx1km[1::2]=0.5*(sgx1km[0:-1:2]+sgx1km[2::2])
sgy1km=np.zeros(2*nj+1);sgy1km[::2]=yq1km;sgy1km[1::2]=0.5*(sgy1km[0:-1:2]+sgy1km[2::2])


print('x bnds: ',(xq1km[0],xq1km[-1]))
print('y bnds: ',(yq1km[0],yq1km[-1]))

fralpha1km=ds['friction_coefficient'].fillna(0.).load().data
Bbar1km=ds['Bbar'].fillna(0.).load().data


grid_1km = quadmesh('Grnld_1km')
print('Greenland supergrid: ',(grid_1km.nj,grid_1km.ni))
grid_2km = quadmesh('Grnld_2km',parent=grid_1km)
print('Greenland 2km supergrid: ',(grid_2km.nj,grid_2km.ni))
grid_4km = quadmesh('Grnld_4km',parent=grid_2km)
print('Greenland 2km supergrid: ',(grid_4km.nj,grid_4km.ni))

grid_1km.add_rheology()
grid_2km.add_rheology(grid_1km)
grid_4km.add_rheology(grid_2km)

for g in [grid_1km,grid_2km,grid_4km]:
    g.print()

grid_1km.to_netcdf()
grid_2km.to_netcdf()
grid_4km.to_netcdf()
raise()

#trans = Transformer.from_crs(bm_proj,m_proj)
#X_bm, Y_bm = trans.transform(Xp_bm,Yp_bm)
