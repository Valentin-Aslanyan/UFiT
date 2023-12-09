

"""
Create a 3D image of Q and some fieldlines 
"""


Q_file="Example.flf"
fieldline_file="Example2.flf"


import sys
UFiT_home_directory='/change/this/path'
sys.path.append(UFiT_home_directory)
from UFiT_Functions_Python import *
from mayavi import mlab


fieldline_flf=read_UFiT_output(UFiTfile_name=fieldline_file)

mlab.figure(bgcolor=(1.0,1.0,1.0),size=(1200,1000))
if fieldline_flf.geometry==0:
	Q_flf=read_UFiT_output(UFiTfile_name=Q_file)
	dim1,dim2,Q=get_Q_2D(Q_flf)
	dim1,dim2=np.meshgrid(dim1,dim2)
	Q[np.isinf(Q)]=1E6
	Q[np.isinf(-Q)]=-1E6
	slogQ=np.sign(Q)*np.log(np.clip(abs(Q),2,1E6))/np.log(10.0)
	dim3=Q_flf.coord3[0]*np.ones(np.shape(dim1))

	for idx in range(len(fieldline_flf.fieldlines)):
		mlab.plot3d(fieldline_flf.fieldlines[idx][:,0],fieldline_flf.fieldlines[idx][:,1],fieldline_flf.fieldlines[idx][:,2],line_width=0.01,color=(0/255,0/255,0/255),tube_radius=0.002)

	plane_mesh=mlab.mesh(dim1,dim2,dim3,scalars=slogQ,colormap='bwr',vmin=-6,vmax=6)

if fieldline_flf.geometry==1:
	max_R=3.0
	Q_flf=read_UFiT_output(UFiTfile_name=Q_file)
	dim1,dim2,Q=get_Q_2D(Q_flf)
	dim1,dim2=np.meshgrid(dim1,dim2)
	Q[np.isinf(Q)]=1E6
	Q[np.isinf(-Q)]=-1E6
	slogQ=np.sign(Q)*np.log(np.clip(abs(Q),2,1E6))/np.log(10.0)
	dim3=Q_flf.coord1[0]*np.ones(np.shape(dim1))

	plane_X=dim3*np.sin(dim1)*np.cos(dim2)
	plane_Y=dim3*np.sin(dim1)*np.sin(dim2)
	plane_Z=dim3*np.cos(dim1)

	for idx in range(len(fieldline_flf.fieldlines)):
		if max(fieldline_flf.fieldlines[idx][:,0])>max_R:
			max_idx=np.argmin(abs(fieldline_flf.fieldlines[idx][:,0]-max_R))
		else:
			max_idx=len(fieldline_flf.fieldlines[idx][:,0])
		fieldline_X=fieldline_flf.fieldlines[idx][:max_idx,0]*np.sin(fieldline_flf.fieldlines[idx][:max_idx,1])*np.cos(fieldline_flf.fieldlines[idx][:max_idx,2])
		fieldline_Y=fieldline_flf.fieldlines[idx][:max_idx,0]*np.sin(fieldline_flf.fieldlines[idx][:max_idx,1])*np.sin(fieldline_flf.fieldlines[idx][:max_idx,2])
		fieldline_Z=fieldline_flf.fieldlines[idx][:max_idx,0]*np.cos(fieldline_flf.fieldlines[idx][:max_idx,1])
		mlab.plot3d(fieldline_X,fieldline_Y,fieldline_Z,line_width=0.01,color=(0/255,0/255,0/255),tube_radius=0.002)

	plane_mesh=mlab.mesh(plane_X,plane_Y,plane_Z,scalars=slogQ,colormap='bwr',vmin=-6,vmax=6)
	mlab.view(azimuth=180.0, elevation=90.0, distance=5.0)

mlab.savefig("Spherical_Example3D.png")
mlab.savefig("Spherical_Example3D.png")
mlab.show()



