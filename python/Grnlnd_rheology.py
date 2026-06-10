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



    def add_rheology(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.fralpha = 0.25*(fralpha1km[:-1,:-1]+fralpha1km[:-1,1:]+fralpha1km[1:,:-1]+fralpha1km[1:,1:])
                self.Bbar = 0.25*(Bbar1km[:-1,:-1]+Bbar1km[:-1,1:]+Bbar1km[1:,:-1]+Bbar1km[1:,1:])
                #self.A_glen = self.Bbar.fillna(168239086.57399723)**(-3)
                #self.tau_b_beta = self.fralpha.fillna(1.e5)**2.0
                Bbar=self.Bbar
                Bbar[Bbar==0.]=168239086.57399723
                self.A_glen = Bbar**(-3)
                self.tau_b_beta = self.fralpha**2.0
        else:
            self.fralpha = 0.25*((parent.fralpha[:-1:2,:-1:2]+parent.fralpha[:-1:2,1::2])+\
                             parent.fralpha[1::2,:-1:2]+parent.fralpha[1::2,1::2])
            print('Added fralpha elevation for child grid mean/std: ',self.fralpha.mean(),self.fralpha.std())
            self.Bbar = 0.25*((parent.Bbar[:-1:2,:-1:2]+parent.Bbar[:-1:2,1::2])+\
                             parent.Bbar[1::2,:-1:2]+parent.Bbar[1::2,1::2])
            print('Added Bbar elevation for child grid mean/std: ',self.Bbar.mean(),self.Bbar.std())
            #self.A_glen = self.Bbar.fillna(168239086.57399723)**(-3)
            #self.tau_b_beta = self.fralpha.fillna(1.e5)**2.0
            self.A_glen = self.Bbar**(-3)
            self.tau_b_beta = self.fralpha**2.0


    def add_vel(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.uvel =  vx_obs1km
                self.vvel =  vy_obs1km
                self.float_frac = float_frac1km
                print('Added uvel for 1km grid mean/std: ',self.uvel.mean(),self.uvel.std())
                print('Added vvel for 1km grid mean/std: ',self.vvel.mean(),self.vvel.std())
        else:
            """
            ## This is the preferred way - mjh
            vh = parent.vx*parent.thickness
            self.vx =0.25*((vh[:-1:2,:-1:2]+vh[:-1:2,1::2])+\
                      vh[1::2,:-1:2]+vh[1::2,1::2])/np.maximum(self.thickness,1.e-10)
            vh = parent.vy*parent.thickness
            self.vy =0.25*((vh[:-1:2,:-1:2]+vh[:-1:2,1::2])+\
                      vh[1::2,:-1:2]+vh[1::2,1::2])/np.maximum(self.thickness,1.e-10)
            """
            self.uvel = parent.uvel[::2,::2]
            self.vvel = parent.vvel[::2,::2]

            print('Added vx for child grid mean/std: ',self.uvel.mean(),self.uvel.std())
            print('Added vy for child grid mean/std: ',self.vvel.mean(),self.vvel.std())

            self.float_frac = 0.25*((parent.float_frac[:-1:2,:-1:2]+parent.float_frac[:-1:2,1::2])+\
                             parent.float_frac[1::2,:-1:2]+parent.float_frac[1::2,1::2])

    def add_velmask(self,parent=None):
        if parent is None:
            if self.level == 1:
                self.umask = umask1km
                print('Added umask for 1km grid max/min: ',self.umask.max(),self.umask.min())
                self.vmask = vmask1km
                print('Added vmask for 1km grid max/min: ',self.vmask.max(),self.vmask.min())
        else:
            self.umask = parent.umask[::2,::2]
            self.vmask = parent.vmask[::2,::2]

            print('Added umask for child grid max/min: ',self.umask.max(),self.umask.min())
            print('Added umask for child grid max/min: ',self.vmask.max(),self.vmask.min())

    def to_netcdf(self,path=None):
        xh=self.x[1::2]
        yh=self.y[1::2]
        xq=self.x[::2]
        yq=self.y[::2]
        path_out = self.name+'_rheology.nc'
        if path is not None: path_out=path
        print('saving to netcdf, ',path_out, self.name)
        dA_fr=xr.DataArray(name='friction_coefficient',data=self.fralpha,dims=["y","x"],coords=dict(y=(["y"],yh),x=(["x"],xh)))
        dA_A_glen=xr.DataArray(name='A_glen',data=self.A_glen,dims=["y","x"],coords=dict(y=(["y"],yh),x=(["x"],xh)))
        dA_bb=xr.DataArray(name='Bbar',data=self.Bbar,dims=["y","x"],coords=dict(y=(["y"],yh),x=(["x"],xh)))
        dA_tau_b_beta=xr.DataArray(name='tau_b_beta',data=self.tau_b_beta,dims=["y","x"],coords=dict(y=(["y"],yh),x=(["x"],xh)))
        for d_ in [dA_fr,dA_A_glen,dA_bb,dA_tau_b_beta]:
            d_.coords['x'].attrs['cartesian_axis']='X'
            d_.coords['y'].attrs['cartesian_axis']='Y'

        ds_out=xr.merge([dA_fr,dA_bb,dA_A_glen,dA_tau_b_beta])
        ds_out.to_netcdf(path_out)
        path_out = self.name+'_velocity.nc'
        dA_uvel=xr.DataArray(name='uvel',data=self.uvel,dims=["yp","xp"],coords=dict(yp=(["yp"],yq),xp=(["xp"],xq)))
        dA_vvel=xr.DataArray(name='vvel',data=self.vvel,dims=["yp","xp"],coords=dict(yp=(["yp"],yq),xp=(["xp"],xq)))
        dA_umask=xr.DataArray(name='umask',data=self.umask,dims=["yp","xp"],coords=dict(yp=(["yp"],yq),xp=(["xp"],xq)))
        dA_vmask=xr.DataArray(name='vmask',data=self.vmask,dims=["yp","xp"],coords=dict(yp=(["yp"],yq),xp=(["xp"],xq)))
        for d_ in [dA_uvel,dA_vvel,dA_umask,dA_vmask]:
            d_.coords['xp'].attrs['cartesian_axis']='X'
            d_.coords['yp'].attrs['cartesian_axis']='Y'

        dA_floatfrac=xr.DataArray(name='float_frac',data=self.float_frac,dims=["y","x"],coords=dict(y=(["y"],yh),x=(["x"],xh)))
        dA_floatfrac.coords['x'].attrs['cartesian_axis']='X'
        dA_floatfrac.coords['y'].attrs['cartesian_axis']='Y'
        ds_out=xr.merge([dA_uvel,dA_vvel,dA_umask,dA_vmask,dA_floatfrac])
        ds_out.to_netcdf(path_out)

#path='INPUT_NS/GreenlandISMIP6_BMa5_Control_drag_weertman_abslog_BMa5_interp_d1.nc'
path='INPUT_NS/GreenlandISMIP7_Control_drag_weertman_abslog_final_interp_d1.nc'
print('Opening dataset: ',path)
ds=xr.open_dataset(path)
#Cell corners xq,yq (m)
xq1km=ds['x'].load().data
yq1km=ds['y'].load().data
nip=xq1km.shape[0];njp=yq1km.shape[0]
ni=nip-1;nj=njp-1;
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

fralpha1km=ds['friction_coefficient'].fillna(0.).load().data
Bbar1km=ds['Bbar'].fillna(0.).load().data
yr_p_sec=1.0/8.64e4/365.25
vx_obs1km=ds['vx_obs'].fillna(0.).load().data*yr_p_sec # convert to m/s from m/yr
vy_obs1km=ds['vy_obs'].fillna(0.).load().data*yr_p_sec
umask1km = np.zeros((njp,nip))-1
vmask1km = np.zeros((njp,nip))-1
float_frac1km = np.zeros((nj,ni))


grid_1km = quadmesh('Grnld_1km')
print('Greenland supergrid: ',(grid_1km.nj,grid_1km.ni))
grid_2km = quadmesh('Grnld_2km',parent=grid_1km)
print('Greenland 2km supergrid: ',(grid_2km.nj,grid_2km.ni))
grid_4km = quadmesh('Grnld_4km',parent=grid_2km)
print('Greenland 2km supergrid: ',(grid_4km.nj,grid_4km.ni))

grid_1km.add_rheology()
grid_2km.add_rheology(grid_1km)
grid_4km.add_rheology(grid_2km)
grid_1km.add_velmask()
grid_2km.add_velmask(grid_1km)
grid_4km.add_velmask(grid_2km)
grid_1km.add_vel()
grid_2km.add_vel(grid_1km)
grid_4km.add_vel(grid_2km)

"""
ds_1km=xr.open_dataset('Grnld_1km.nc')
ds_2km=xr.open_dataset('Grnld_2km.nc')
ds_4km=xr.open_dataset('Grnld_4km.nc')

mask_1km=ds_1km['mask'].load().data
thick_1km=ds_1km['thickness'].load().data
mask_1km[thick_1km>0.]=1.0
ds_1km['mask'].data=mask_1km
mask_2km=ds_2km['mask'].load().data
thick_2km=ds_2km['thickness'].load().data
mask_2km[thick_2km>0.]=1.0
ds_2km['mask'].data=mask_2km
mask_4km=ds_4km['mask'].load().data
thick_4km=ds_4km['thickness'].load().data
mask_4km[thick_4km>0.]=1.0
ds_4km['mask'].data=mask_4km

ds_1km_=ds_1km.rename({'mask':'h_mask'})
ds_2km_=ds_2km.rename({'mask':'h_mask'})
ds_4km_=ds_4km.rename({'mask':'h_mask'})

ds_1km_.to_netcdf('Grnld_hmask_1km.nc',mode='w')
ds_2km_.to_netcdf('Grnld_hmask_2km.nc',mode='w')
ds_4km_.to_netcdf('Grnld_hmask_4km.nc',mode='w')
"""
for g in [grid_1km,grid_2km,grid_4km]:
    g.print()

grid_1km.to_netcdf()
grid_2km.to_netcdf()
grid_4km.to_netcdf()



raise()

#trans = Transformer.from_crs(bm_proj,m_proj)
#X_bm, Y_bm = trans.transform(Xp_bm,Yp_bm)
