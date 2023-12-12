

"""
Download a synoptic HMI magnetogram and calculate the squashing factor Q by calling UFiT directly from Python
"""

target_CR=2112	#Carrington rotation

num_r=60	#dimensions of B grid
num_t=180
num_p=360

Rss=2.5		#Source surface

numin_r=1	#Dimensions of output
numin_t=600
numin_p=1200



import sys
UFiT_home_directory='/change/this/path'
sys.path.append(UFiT_home_directory)
from UFiT_Functions_Python import *
compiled_fortran_path=os.path.join(UFiT_home_directory,"UFiT_Python_Callable.so")
import matplotlib.pyplot as plt
import wget
from astropy.io import fits
import astropy.units
import pfsspy
import sunpy


#limits of output grid
r_in_limits=[1.05,1.05]
theta_in_limits=[0.0,np.pi]
phi_in_limits=[0.0,2.0*np.pi]


def get_HMI_magnetogram(cr_target,data_dir):
	if cr_target<2096:
		print("CR "+str(cr_target)+" too early for HMI")
		return ""
	target_filename='hmi.Synoptic_Mr.'+str(cr_target)+'.fits'
	
	if target_filename not in os.listdir(data_dir):
		wget.download('http://jsoc.stanford.edu/data/hmi/synoptic/'+target_filename,out=data_dir)
	if target_filename in os.listdir(data_dir):
		return os.path.join(data_dir,target_filename)
	else:
		print("Could not get HMI map for CR "+str(cr_target))
		return ""


#Download magnetogram fits file
HMI_filename=get_HMI_magnetogram(target_CR,UFiT_home_directory)
HMI_map = sunpy.map.Map(HMI_filename)

#Downsample and remove NaNs as required by pfsspy
HMI_map = HMI_map.resample([num_p, num_t] * astropy.units.pix)
HMI_map.data[np.isnan(HMI_map.data)]=0.0	#NaNs set to zero

#Compute B from PFSS and transform the field and coordinates into required form
pfss_in = pfsspy.Input(HMI_map, num_r, Rss)
pfss_out = pfsspy.pfss(pfss_in)
grid_r=np.exp(pfss_out.grid.rg)
grid_t=np.arccos(pfss_out.grid.sg)
grid_p=pfss_out.grid.pg
B=np.zeros((3,len(grid_r),len(grid_t),len(grid_p)))
B[0,:,:,:]=pfss_out.bg[:,:,:,2].transpose(2,1,0)
B[1,:,:,:]=pfss_out.bg[:,:,:,1].transpose(2,1,0)
B[2,:,:,:]=pfss_out.bg[:,:,:,0].transpose(2,1,0)

#Create fieldline starts
r_in=np.linspace(r_in_limits[0],r_in_limits[1],num=numin_r)
theta_in=np.linspace(theta_in_limits[0],theta_in_limits[1],num=numin_t)
phi_in=np.linspace(phi_in_limits[0],phi_in_limits[1],num=numin_p)

#Create and populate an object with required information; this corresponds to the command line options, input file and B file
UFiT_run=UFiT_call_input()
UFiT_run.geometry=1
UFiT_run.input_type=1
UFiT_run.num_proc=4
UFiT_run.periodic_PHI=True
UFiT_run.normalized_B=True
UFiT_run.include_curvature=False
UFiT_run.grid_regular=True
UFiT_run.load_B=False
UFiT_run.save_fieldlines=False
UFiT_run.save_Q=True
UFiT_run.save_endpoints=False
UFiT_run.write_output = False
UFiT_run.return_output = True
UFiT_run.grid1=grid_r
UFiT_run.grid2=grid_t
UFiT_run.grid3=grid_p
UFiT_run.B_grid=B
UFiT_run.coord1=r_in
UFiT_run.coord2=theta_in
UFiT_run.coord3=phi_in

#Run UFiT
UFiT_run=call_UFiT(compiled_fortran_path,UFiT_run)

#Transform raw output and display
dim1,dim2,Q=get_Q_2D(UFiT_run)
Q[np.isinf(Q)]=1E6
Q[np.isinf(-Q)]=-1E6
slogQ=np.sign(Q)*np.log(np.clip(abs(Q),2,1E6))/np.log(10.0)

plt.figure(figsize=(10,6))
if UFiT_run.geometry==0:
	plt.pcolormesh(dim1,dim2,slogQ,cmap="bwr",vmin=-6,vmax=6,rasterized=True)
if UFiT_run.geometry==1:
	plt.pcolormesh(dim2*180.0/np.pi,90.0-dim1*180.0/np.pi,np.transpose(slogQ),cmap="bwr",vmin=-6,vmax=6,rasterized=True)
cbar=plt.colorbar()
cbar.solids.set_rasterized(True)
cbar.ax.tick_params(labelsize=20,direction='in', bottom=True, top=True, left=True, right=True)
cbar.set_label(label="$\mathrm{slog}(Q)$",fontsize=20)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
plt.tick_params(axis='both', which='major',labelsize=20,direction='in',bottom=True, top=True, left=True, right=True)
plt.savefig("HMI_Example.png", dpi=100,bbox_inches='tight',pad_inches=0.1)
plt.show()











