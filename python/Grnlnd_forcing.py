import numpy as np  # numerical library
import xarray as xr  # netCDF library
from pyproj import CRS
from pyproj.enums import WktVersion
from pyproj import Transformer


class quadmesh:
    """
    Rectangular grid class
    """
    def __init__(self,name=None,parent=None):
        self.name = name
        if parent==None:
            self.x = sgx1km
            self.y = sgy1km
            self.level = 1
        else:
            self.x = parent.x[::2]
            self.y = parent.y[::2]
            self.level = parent.level + 1
        self.ni=int(0.5*(self.x.shape[0]-1));self.nj=int(0.5*(self.y.shape[0]-1))

    def print(self):
        print(self.name,' grid dimensions: ',(self.nj,self.ni))
        try:
            print('friction_coefficient : (max/min)',self.fralpha.shape,self.fralpha.max(),self.fralpha.min())
        except:
            pass
        try:
            print('Bbar : (max/min)',self.Bbar.shape,self.Bbar.max(),self.Bbar.min())
        except:
            pass



    def add_forcing(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.acabf = 0.25*(acabf1km[:,:-1,:-1]+acabf1km[:,:-1,1:]+acabf1km[:,1:,:-1]+acabf1km[:,1:,1:])
            print('Added acabf for parent grid mean/std: ',self.acabf.mean(),self.acabf.std())

        else:
            self.acabf = 0.25*((parent.acabf[:,:-1:2,:-1:2]+parent.acabf[:,:-1:2,1::2])+\
                             parent.acabf[:,1::2,:-1:2]+parent.acabf[:,1::2,1::2])
            print('Added acabf for child grid mean/std: ',self.acabf.mean(),self.acabf.std())


    def to_netcdf(self,path=None):
        xh=self.x[1::2]
        yh=self.y[1::2]
        path_out = self.name+'_forcing.nc'
        if path is not None: path_out=path
        print('saving to netcdf, ',path_out, self.name)
        dA_acabf=xr.DataArray(name='acabf',data=self.acabf,dims=["time","y","x"],coords=dict(time=(["time"],time.data),y=(["y"],yh),x=(["x"],xh)))
        dA_acabf.time.attrs['units']=time.units
        dA_acabf.time.attrs['calendar']='noleap'
        dA_acabf.time.attrs['axis']='T'
        dA_acabf.time.attrs['modulo']=' '
        dA_acabf.time.attrs['cartesian_axis']='T'
        dA_acabf.x.attrs['cartesian_axis']='X'
        dA_acabf.y.attrs['cartesian_axis']='Y'
        dA_acabf.to_netcdf(path_out,unlimited_dims=['time',])

path='FORCING/acabf_GrIS_CESM2-WACCM_historical_SDBN1_v2_1850.nc'
print('Opening dataset: ',path)
ds=xr.open_dataset(path,decode_times=False)
#Cell corners xq,yq (m)
xq1km=ds['x'].load().data
yq1km=ds['y'].load().data
#print('flipping y-axis for MOM6')
#yq1km=yq1km[::-1]

ni=xq1km.shape[0]-1;nj=yq1km.shape[0]-1
print('Grid size= ',(nj,ni))
print('x-node coord range: ',(xq1km[0],xq1km[-1]))
print('y-node coord range: ',(yq1km[0],yq1km[-1]))
xh1km = 0.5*(xq1km[:-1]+xq1km[1:])
yh1km = 0.5*(yq1km[:-1]+yq1km[1:])
#proj4text=ds['mapping'].proj4text
#proj = CRS.from_proj4(proj4text)
#print(proj.to_wkt(WktVersion.WKT1_GDAL, pretty=True))
# Distance between cell centers dx,dy (m)
dx=xq1km[1:]-xq1km[:-1]
dy=yq1km[1:]-yq1km[:-1]
xbnds = (xq1km[0],xq1km[-1])
ybnds = (yq1km[0],yq1km[-1])
# Supergrid (h + q)
sgx1km=np.zeros(2*ni+1);sgx1km[::2]=xq1km;sgx1km[1::2]=0.5*(sgx1km[0:-1:2]+sgx1km[2::2])
sgy1km=np.zeros(2*nj+1);sgy1km[::2]=yq1km;sgy1km[1::2]=0.5*(sgy1km[0:-1:2]+sgy1km[2::2])


print('x bnds: ',(xq1km[0],xq1km[-1]))
print('y bnds: ',(yq1km[0],yq1km[-1]))

acabf1km=ds['acabf'].fillna(0.).load().data
#flip y axis
acabf1km=acabf1km[::-1,:]

time = ds['time']


grid_1km = quadmesh('Grnld_1km')
print('Greenland supergrid: ',(grid_1km.nj,grid_1km.ni))
grid_2km = quadmesh('Grnld_2km',parent=grid_1km)
print('Greenland 2km supergrid: ',(grid_2km.nj,grid_2km.ni))
grid_4km = quadmesh('Grnld_4km',parent=grid_2km)
print('Greenland 2km supergrid: ',(grid_4km.nj,grid_4km.ni))

grid_1km.add_forcing()
grid_2km.add_forcing(grid_1km)
grid_4km.add_forcing(grid_2km)

for g in [grid_1km,grid_2km,grid_4km]:
    g.print()

grid_1km.to_netcdf()
grid_2km.to_netcdf()
grid_4km.to_netcdf()
raise()

#trans = Transformer.from_crs(bm_proj,m_proj)
#X_bm, Y_bm = trans.transform(Xp_bm,Yp_bm)
