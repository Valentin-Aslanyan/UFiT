

visualize_3D = False
l_scale = 1.0
r_max = 3.0
z_lims = [-2.0,2.0]
num_r = 200
num_p = 120
num_z = 100
B_z_0 = 1.0
j0 = 1.0


view_pitch = 40.0
view_azimuth = 0.0
view_elevation = 0.0
view_roll = 0.0
view_distance = 35.0
view_focalpoint = (r_max,0.0,z_lims[0])


import sys
UFiT_home_directory='/change/this/path'
sys.path.append(UFiT_home_directory)
from UFiT_Functions_Python import *
from USlip_Functions_Python import *
compiled_UFiT_path=os.path.join(UFiT_home_directory,"UFiT_Python_Callable.so")
compiled_USlip_path=os.path.join(UFiT_home_directory,"USlip_Python_Callable.so")
import matplotlib.pyplot as plt


####//\\  Set up analytic grids  //\\

#Set r grid to start just ABOVE 0, or there is a numerical singularity (not analytically, though)
r_cyl_grid=np.linspace(0.01,r_max,num=num_r)
p_cyl_grid=np.linspace(0.0,2.0*np.pi*(num_p-1)/num_p,num=num_p)
z_cyl_grid=np.linspace(z_lims[0],z_lims[1],num=num_z)

x_cart_grid=np.linspace(-r_max,r_max,num=num_r)
y_cart_grid=np.linspace(-r_max,r_max,num=num_r)
z_cart_grid=np.linspace(z_lims[0],z_lims[1],num=num_z)

B_p_C = 0.5*j0*l_scale**2
B_cyl = np.zeros((3,num_r,num_p,num_z))
B_cyl[1,:,:,:] = (B_p_C*(1.0-np.exp(-(r_cyl_grid/l_scale)**2))/r_cyl_grid)[:,None,None]*np.ones((num_p))[None,:,None]*np.ones((num_z))[None,None,:]
B_cyl[2,:,:,:] = B_z_0

j_cyl_analytic = np.zeros((3,num_r,num_p,num_z))
j_cyl_analytic[2,:,:,:] = (B_p_C*2.0/l_scale**2*np.exp(-(r_cyl_grid/l_scale)**2))[:,None,None]*np.ones((num_p))[None,:,None]*np.ones((num_z))[None,None,:]

Sigma_cyl_analytic = np.zeros((3,num_r,num_p,num_z))
W_cyl=2.0*j0/l_scale**2*(r_cyl_grid*np.exp(-(r_cyl_grid/l_scale)**2))[:,None,None]*np.ones((num_p))[None,:,None]*np.ones((num_z))[None,None,:]
Sigma_cyl_analytic[1,:,:,:] = -W_cyl*B_z_0**2/(B_cyl[1,:,:,:]**2+B_z_0**2)**1.5
Sigma_cyl_analytic[2,:,:,:] = W_cyl*B_z_0*B_cyl[1,:,:,:]/(B_cyl[1,:,:,:]**2+B_z_0**2)**1.5


B_cart=np.zeros((3,num_r,num_r,num_z))
xmg_cart,ymg_cart=np.meshgrid(x_cart_grid,y_cart_grid,indexing='ij')
B_cart[0,:,:,:] = -(B_p_C*(1.0-np.exp(-(xmg_cart**2+ymg_cart**2)/l_scale**2))/(xmg_cart**2+ymg_cart**2)*ymg_cart)[:,:,None]*np.ones((num_z))[None,None,:]
B_cart[1,:,:,:] = (B_p_C*(1.0-np.exp(-(xmg_cart**2+ymg_cart**2)/l_scale**2))/(xmg_cart**2+ymg_cart**2)*xmg_cart)[:,:,None]*np.ones((num_z))[None,None,:]
B_cart[2,:,:,:] = B_z_0

j_cart_analytic = np.zeros((3,num_r,num_r,num_z))
j_cart_analytic[2,:,:,:] = (B_p_C*2.0/l_scale**2*np.exp(-(xmg_cart**2+ymg_cart**2)/l_scale**2))[:,:,None]*np.ones((num_z))[None,None,:]

Sigma_cart_analytic = np.zeros((3,num_r,num_r,num_z))
W_cart=2.0*j0/l_scale**2*(np.sqrt(xmg_cart**2+ymg_cart**2)*np.exp(-(xmg_cart**2+ymg_cart**2)/l_scale**2))[:,:,None]*np.ones((num_z))[None,None,:]
Sigma_cart_analytic[0,:,:,:] = (W_cart*B_z_0**2/(B_cart[0,:,:,:]**2+B_cart[1,:,:,:]**2+B_z_0**2)**1.5)*((ymg_cart/np.sqrt(xmg_cart**2+ymg_cart**2))[:,:,None]*np.ones((num_z))[None,None,:])
Sigma_cart_analytic[1,:,:,:] = -(W_cart*B_z_0**2/(B_cart[0,:,:,:]**2+B_cart[1,:,:,:]**2+B_z_0**2)**1.5)*((xmg_cart/np.sqrt(xmg_cart**2+ymg_cart**2))[:,:,None]*np.ones((num_z))[None,None,:])
Sigma_cart_analytic[2,:,:,:] = W_cart*B_z_0*np.sqrt(B_cart[0,:,:,:]**2+B_cart[1,:,:,:]**2)/(B_cart[0,:,:,:]**2+B_cart[1,:,:,:]**2+B_z_0**2)**1.5

####//\\                         //\\


####|/\|       Call USlip        |/\|

USlip_run=USlip_call_input()
USlip_run.geometry=2
USlip_run.num_proc=4
USlip_run.periodic_PHI=True
USlip_run.grid_regular=True
USlip_run.load_B=False
USlip_run.write_output = False
USlip_run.return_output = True
USlip_run.grid1=r_cyl_grid
USlip_run.grid2=p_cyl_grid
USlip_run.grid3=z_cyl_grid
USlip_run.B_grid=B_cyl
USlip_run=call_USlip(compiled_USlip_path,USlip_run)

USlip_run2=USlip_call_input()
USlip_run2.geometry=0
USlip_run2.num_proc=4
USlip_run2.periodic_PHI=True
USlip_run2.grid_regular=True
USlip_run2.load_B=False
USlip_run2.write_output = False
USlip_run2.return_output = True
USlip_run2.grid1=x_cart_grid
USlip_run2.grid2=y_cart_grid
USlip_run2.grid3=z_cart_grid
USlip_run2.B_grid=B_cart
USlip_run2=call_USlip(compiled_USlip_path,USlip_run2)

####|/\|                         |/\|


## 2D plots with matplotlib
xmg_cyl = np.outer(r_cyl_grid,np.cos(p_cyl_grid))
ymg_cyl = np.outer(r_cyl_grid,np.sin(p_cyl_grid))

fig = plt.figure(figsize=(30,10))

plt.subplot(2,6,1)
plt.pcolormesh(xmg_cyl,ymg_cyl,B_cyl[1,:,:,0]**2)
plt.title("Cylindrical, $B_\\phi^2$")

plt.subplot(2,6,2)
plt.pcolormesh(xmg_cyl,ymg_cyl,j_cyl_analytic[2,:,:,0])
plt.title("Cylindrical, $j_z$, analytic")

plt.subplot(2,6,3)
plt.pcolormesh(xmg_cyl,ymg_cyl,USlip_run.j_grid[2,:,:,0])
plt.title("Cylindrical, $j_z$, from USlip")

plt.subplot(2,6,4)
plt.pcolormesh(xmg_cyl,ymg_cyl,Sigma_cyl_analytic[0,:,:,0]**2+Sigma_cyl_analytic[1,:,:,0]**2+Sigma_cyl_analytic[2,:,:,0]**2)
plt.title("Cylindrical, $|\\Sigma|^2$, analytic")

plt.subplot(2,6,5)
plt.pcolormesh(xmg_cyl,ymg_cyl,USlip_run.sigma_grid[0,:,:,0]**2+USlip_run.sigma_grid[1,:,:,0]**2+USlip_run.sigma_grid[2,:,:,0]**2)
plt.title("Cylindrical, $|\\Sigma|^2$, from Uslip")

plt.subplot(2,6,6)
plt.pcolormesh(xmg_cyl,ymg_cyl,USlip_run.sigmaalpha_grid[0,:,:,0]**2+USlip_run.sigmaalpha_grid[1,:,:,0]**2+USlip_run.sigmaalpha_grid[2,:,:,0]**2)
plt.title("Cylindrical, $|\\Sigma_{\\alpha}|^2$")

plt.subplot(2,6,7)
plt.pcolormesh(xmg_cart,ymg_cart,B_cart[0,:,:,0]**2+B_cart[1,:,:,0]**2)
plt.title("Cartesian, $B_x^2+B_y^2\\equiv B_\\phi^2$")

plt.subplot(2,6,8)
plt.pcolormesh(xmg_cart,ymg_cart,j_cart_analytic[2,:,:,0])
plt.title("Cartesian, $j_z$, analytic")

plt.subplot(2,6,9)
plt.pcolormesh(xmg_cart,ymg_cart,USlip_run2.j_grid[2,:,:,0])
plt.title("Cartesian, $j_z$, from USlip")

plt.subplot(2,6,10)
plt.pcolormesh(xmg_cart,ymg_cart,Sigma_cart_analytic[0,:,:,0]**2+Sigma_cart_analytic[1,:,:,0]**2+Sigma_cart_analytic[2,:,:,0]**2)
plt.title("Cartesian, $|\\Sigma|^2$, analytic")

plt.subplot(2,6,11)
plt.pcolormesh(xmg_cart,ymg_cart,USlip_run2.sigma_grid[0,:,:,0]**2+USlip_run2.sigma_grid[1,:,:,0]**2+USlip_run2.sigma_grid[2,:,:,0]**2)
plt.title("Cartesian, $|\\Sigma|^2$, from Uslip")

plt.subplot(2,6,12)
plt.pcolormesh(xmg_cart,ymg_cart,USlip_run2.sigmaalpha_grid[0,:,:,0]**2+USlip_run2.sigmaalpha_grid[1,:,:,0]**2+USlip_run2.sigmaalpha_grid[2,:,:,0]**2)
plt.title("Cartesian, $|\\Sigma_{\\alpha}|^2$")

plt.savefig("Gaussian_Current_USlip.png", dpi=100,bbox_inches='tight',pad_inches=0.1)


## 3D plots with mayavi
if visualize_3D:
	from mayavi import mlab

	fl_per_circle=8
	fl_starts_r=np.concatenate((0.3*r_max*np.ones((fl_per_circle)),0.6*r_max*np.ones((fl_per_circle))))
	fl_starts_p=np.concatenate((np.linspace(0,2.0*np.pi*(fl_per_circle-1)/fl_per_circle,num=fl_per_circle),np.linspace(0,2.0*np.pi*(fl_per_circle-1)/fl_per_circle,num=fl_per_circle)))
	fl_starts_z=z_lims[0]*np.ones((len(fl_starts_r)))

	fl_starts_x=fl_starts_r*np.cos(fl_starts_p)
	fl_starts_y=fl_starts_r*np.sin(fl_starts_p)

	UFiT_run=UFiT_call_input()
	UFiT_run.geometry=2
	UFiT_run.input_type=0
	UFiT_run.num_proc=4
	UFiT_run.periodic_PHI=True
	UFiT_run.normalized_B=True
	UFiT_run.include_curvature=True
	UFiT_run.grid_regular=True
	UFiT_run.load_B=False
	UFiT_run.save_fieldlines=True
	UFiT_run.save_Q=False
	UFiT_run.save_endpoints=False
	UFiT_run.write_output = False
	UFiT_run.return_output = True
	UFiT_run.coord1=fl_starts_r
	UFiT_run.coord2=fl_starts_p
	UFiT_run.coord3=fl_starts_z
	UFiT_run.grid1=r_cyl_grid
	UFiT_run.grid2=p_cyl_grid
	UFiT_run.grid3=z_cyl_grid
	UFiT_run.B_grid=B_cyl
	UFiT_run=call_UFiT(compiled_UFiT_path,UFiT_run)

	UFiT_run2=UFiT_call_input()
	UFiT_run2.geometry=0
	UFiT_run2.input_type=0
	UFiT_run2.num_proc=4
	UFiT_run2.normalized_B=True
	UFiT_run2.include_curvature=True
	UFiT_run2.grid_regular=True
	UFiT_run2.load_B=False
	UFiT_run2.save_fieldlines=True
	UFiT_run2.save_Q=False
	UFiT_run2.save_endpoints=False
	UFiT_run2.write_output = False
	UFiT_run2.return_output = True
	UFiT_run2.coord1=fl_starts_x
	UFiT_run2.coord2=fl_starts_y
	UFiT_run2.coord3=fl_starts_z
	UFiT_run2.grid1=x_cart_grid
	UFiT_run2.grid2=y_cart_grid
	UFiT_run2.grid3=z_cart_grid
	UFiT_run2.B_grid=B_cart
	UFiT_run2=call_UFiT(compiled_UFiT_path,UFiT_run2)


	fieldlines_cyl2cart=[]
	fieldlines_cart2cart=[]
	for idx in range(np.shape(UFiT_run.fieldline_allpos)[0]):
		fieldlines_cyl2cart.append(UFiT_run.fieldline_allpos[idx,UFiT_run.fieldline_pts[idx]-1:UFiT_run.fieldline_pts[idx]+UFiT_run.fieldline_ptn[idx]-1,:])
	for idx in range(np.shape(UFiT_run2.fieldline_allpos)[0]):
		fieldlines_cart2cart.append(UFiT_run2.fieldline_allpos[idx,UFiT_run2.fieldline_pts[idx]-1:UFiT_run2.fieldline_pts[idx]+UFiT_run2.fieldline_ptn[idx]-1,:])


	fig=mlab.figure(bgcolor=(1.0,1.0,1.0),size=(1200,1000))

	mlab.clf()

	mlab.view(azimuth=view_azimuth, elevation=view_elevation, roll=view_roll, distance=view_distance, focalpoint=view_focalpoint)
	camera = fig.scene.camera

	plane_mesh=mlab.mesh(xmg_cyl,ymg_cyl,0.0*xmg_cyl+z_lims[0],scalars=USlip_run.sigma_grid[0,:,:,0]**2+USlip_run.sigma_grid[1,:,:,0]**2+USlip_run.sigma_grid[2,:,:,0]**2)#,colormap='bwr',vmin=-6,vmax=6)
	plane_mesh2=mlab.mesh(xmg_cart+2.1*r_max,ymg_cart,0.0*xmg_cart+z_lims[0],scalars=USlip_run2.sigma_grid[0,:,:,0]**2+USlip_run2.sigma_grid[1,:,:,0]**2+USlip_run2.sigma_grid[2,:,:,0]**2)#,colormap='bwr',vmin=-6,vmax=6)

	mlab.savefig("Gaussian_Current_USlip3D.png")

	mlab.clf()

	mlab.view(azimuth=view_azimuth, elevation=view_elevation, roll=view_roll, distance=view_distance, focalpoint=view_focalpoint)
	camera = fig.scene.camera

	plane_mesh=mlab.mesh(xmg_cyl,ymg_cyl,0.0*xmg_cyl+z_lims[0],scalars=USlip_run.sigma_grid[0,:,:,0]**2+USlip_run.sigma_grid[1,:,:,0]**2+USlip_run.sigma_grid[2,:,:,0]**2)#,colormap='bwr',vmin=-6,vmax=6)
	plane_mesh2=mlab.mesh(xmg_cart+2.1*r_max,ymg_cart,0.0*xmg_cart+z_lims[0],scalars=USlip_run2.sigma_grid[0,:,:,0]**2+USlip_run2.sigma_grid[1,:,:,0]**2+USlip_run2.sigma_grid[2,:,:,0]**2)#,colormap='bwr',vmin=-6,vmax=6)

	mlab.savefig("Gaussian_Current_USlip3D.png")

	mlab.clf()


	for idx in range(len(fieldlines_cart2cart)):
		fieldline_X=fieldlines_cart2cart[idx][:,0]+2.1*r_max
		fieldline_Y=fieldlines_cart2cart[idx][:,1]
		fieldline_Z=fieldlines_cart2cart[idx][:,2]
		mlab.plot3d(fieldline_X,fieldline_Y,fieldline_Z,line_width=0.01,color=(0/255,0/255,0/255),tube_radius=0.02)

	for idx in range(len(fieldlines_cyl2cart)):
		fieldline_X=fieldlines_cyl2cart[idx][:,0]*np.cos(fieldlines_cyl2cart[idx][:,1])
		fieldline_Y=fieldlines_cyl2cart[idx][:,0]*np.sin(fieldlines_cyl2cart[idx][:,1])
		fieldline_Z=fieldlines_cyl2cart[idx][:,2]
		mlab.plot3d(fieldline_X,fieldline_Y,fieldline_Z,line_width=0.01,color=(0/255,0/255,0/255),tube_radius=0.02)


	plane_mesh=mlab.mesh(xmg_cyl,ymg_cyl,0.0*xmg_cyl+z_lims[0],scalars=USlip_run.sigma_grid[0,:,:,0]**2+USlip_run.sigma_grid[1,:,:,0]**2+USlip_run.sigma_grid[2,:,:,0]**2)#,colormap='bwr',vmin=-6,vmax=6)
	plane_mesh2=mlab.mesh(xmg_cart+2.1*r_max,ymg_cart,0.0*xmg_cart+z_lims[0],scalars=USlip_run2.sigma_grid[0,:,:,0]**2+USlip_run2.sigma_grid[1,:,:,0]**2+USlip_run2.sigma_grid[2,:,:,0]**2)#,colormap='bwr',vmin=-6,vmax=6)
	camera = fig.scene.camera

	fig.scene.light_manager.number_of_lights = 5
	fig.scene.light_manager.lights[0].azimuth=0.0
	fig.scene.light_manager.lights[0].elevation=30.0
	fig.scene.light_manager.lights[0].intensity=0.8
	fig.scene.light_manager.lights[0].activate=True
	fig.scene.light_manager.lights[1].azimuth=0.0
	fig.scene.light_manager.lights[1].elevation=70.0
	fig.scene.light_manager.lights[1].intensity=0.4
	fig.scene.light_manager.lights[1].activate=True
	fig.scene.light_manager.lights[2].azimuth=0.0
	fig.scene.light_manager.lights[2].elevation=-70.0
	fig.scene.light_manager.lights[2].intensity=0.3
	fig.scene.light_manager.lights[2].activate=True
	fig.scene.light_manager.lights[3].azimuth=-70.0
	fig.scene.light_manager.lights[3].elevation=0.0
	fig.scene.light_manager.lights[3].intensity=0.3
	fig.scene.light_manager.lights[3].activate=True
	fig.scene.light_manager.lights[4].azimuth=70.0
	fig.scene.light_manager.lights[4].elevation=0.0
	fig.scene.light_manager.lights[4].intensity=0.2
	fig.scene.light_manager.lights[4].activate=True

	mlab.view(azimuth=view_azimuth, elevation=view_elevation, roll=view_roll, distance=view_distance, focalpoint=view_focalpoint)
	mlab.pitch(view_pitch)
	mlab.view(focalpoint=view_focalpoint)


	mlab.savefig("Gaussian_Current_USlip3D.png")


	mlab.show()





