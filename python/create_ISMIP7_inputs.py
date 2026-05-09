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
            print('bed elevation : (max/min/mean)',self.bed.shape,self.bed.max(),self.bed.min(),self.bed.mean())
        except:
            pass
        try:
            print('thickness : (max/min/mean)',self.thickness.shape,self.thickness.max(),self.thickness.min(),self.thickness.mean())
        except:
            pass
        try:
            print('geiod : (max/min/mean)',self.geiod.shape,self.geiod.max(),self.geiod.min(),self.geiod.mean())
        except:
            pass
        try:
            print('mask : (max/min)',self.mask.shape,self.mask.max(),self.mask.min())
        except:
            pass



    def add_bed(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.bed = bed1km[:-1,:-1] # cropping the odd cells off the upper-right side of the array
                print('Added Bed elevation for 1km grid mean/std: ',self.bed.mean(),self.bed.std())
        else:
            self.bed = 0.25*((parent.bed[:-1:2,:-1:2]+parent.bed[:-1:2,1::2])+\
                             parent.bed[1::2,:-1:2]+parent.bed[1::2,1::2])
            print('Added bed elevation for child grid mean/std: ',self.bed.mean(),self.bed.std())

    def add_thickness(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.thickness = thick1km[:-1,:-1] # cropping the odd cells off the upper-right side of the array
                print('Added shelf thickness for 1km grid mean/std: ',self.thickness.mean(),self.thickness.std())
        else:
            self.thickness = 0.25*((parent.thickness[:-1:2,:-1:2]+parent.thickness[:-1:2,1::2])+\
                             parent.thickness[1::2,:-1:2]+parent.thickness[1::2,1::2])
            print('Added thickness for child grid mean/std: ',self.thickness.mean(),self.thickness.std())

    def add_geoid(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.geoid = geoid1km[:-1,:-1] # cropping the odd cells off the upper-right side of the array
                print('Added Geoid for 1km grid mean/std: ',self.geoid.mean(),self.geoid.std())
        else:
            self.geoid = 0.25*((parent.geoid[:-1:2,:-1:2]+parent.geoid[:-1:2,1::2])+\
                             parent.geoid[1::2,:-1:2]+parent.geoid[1::2,1::2])
            print('Added geoid for child grid mean/std: ',self.geoid.mean(),self.geoid.std())

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
        path_out = self.name+'.nc'
        if path is not None: path_out=path
        ds_out = xr.Dataset(data_vars=dict(bed=(["y","x"],self.bed),thickness=(["y","x"],self.thickness),geoid=(["y","x"],self.geoid),mask=(["y","x"],self.mask)),coords=dict(y=yh,x=xh,),attrs=dict(name=self.name),)
        ds_out.to_netcdf(path_out)

path='../INPUT/GreenlandObsISMIP7-v1.3.nc'
print('Opening dataset: ',path)
ds=xr.open_dataset(path)
#Cell centers xh,yh (m)
xh=ds['x1km'].load().data
yh=ds['y1km'].load().data
ni=xh.shape[0];nj=yh.shape[0]
print('Grid size= ',(ni,nj))
print('x center coord range: ',(xh[0],xh[-1]))
print('y center coord range: ',(yh[0],yh[-1]))
proj4text=ds['mapping'].proj4text
proj = CRS.from_proj4(proj4text)
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


# bed depth (m)
xb=ds['x'].load().data
yb=ds['y'].load().data
print('xbed bnds: ',(xb[0],xb[-1]))
print('ybed bnds: ',(yb[0],yb[-1]))
# Find start and end indices of bed data within the 1km domain
[(istart,jstart),(iend,jend)] = find_in_domain([(xb[0],yb[0]),(xb[-1],yb[-1])],xq1km,yq1km)
print('Bed data x-index range: ',(istart,iend))
print('Bed data y-index range: ',(jstart,jend))

bed1km=np.zeros((nj,ni))
thick1km=np.zeros((nj,ni))
geoid1km=np.zeros((nj,ni))
mask1km=np.zeros((nj,ni)).astype('uint8')

bed=ds['bed'].load().data
thick=ds['thickness'].load().data
geoid=ds['geoid'].load().data
mask=ds['mask'].load().data


for j in  np.arange(jstart,jend-1):
    for i in np.arange(istart,iend-1):
        [(i0,j0),(i1,j1)]= find_in_domain([(xq1km[i],yq1km[j]),(xq1km[i+1],yq1km[j+1])],xb,yb)
        bed1km[j,i]=bed[j0:j1,i0:i1].mean()
        thick1km[j,i]=thick[j0:j1,i0:i1].mean()
        geoid1km[j,i]=geoid[j0:j1,i0:i1].mean()
        mask1km[j,i]=np.round(mask[j0:j1,i0:i1].mean()).astype('uint8')

grid_1km = quadmesh('Grnld_1km')
print('Greenland supergrid: ',(grid_1km.nj,grid_1km.ni))
grid_2km = quadmesh('Grnld_2km',parent=grid_1km)
print('Greenland 2km supergrid: ',(grid_2km.nj,grid_2km.ni))
grid_4km = quadmesh('Grnld_4km',parent=grid_2km)
print('Greenland 2km supergrid: ',(grid_4km.nj,grid_4km.ni))

grid_1km.add_bed()
grid_2km.add_bed(grid_1km)
grid_4km.add_bed(grid_2km)
grid_1km.add_thickness()
grid_2km.add_thickness(grid_1km)
grid_4km.add_thickness(grid_2km)
grid_1km.add_geoid()
grid_2km.add_geoid(grid_1km)
grid_4km.add_geoid(grid_2km)
grid_1km.add_mask()
grid_2km.add_mask(grid_1km)
grid_4km.add_mask(grid_2km)

for g in [grid_1km,grid_2km,grid_4km]:
    g.print()

grid_1km.to_netcdf()
grid_2km.to_netcdf()
grid_4km.to_netcdf()
raise()

#trans = Transformer.from_crs(bm_proj,m_proj)
#X_bm, Y_bm = trans.transform(Xp_bm,Yp_bm)
