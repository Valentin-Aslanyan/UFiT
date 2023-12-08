

fieldline_file='Example2.flf'


import sys
UFiT_home_directory='/change/this/path'
sys.path.append(UFiT_home_directory)
from UFiT_Functions_Python import *
from mayavi import mlab


fieldline_flf=read_UFiT_output(UFiTfile_name=fieldline_file)

if fieldline_flf.geometry==0:
	for idx in range(len(fieldline_flf.fieldlines)):
		mlab.plot3d(fieldline_flf.fieldlines[idx][:,0],fieldline_flf.fieldlines[idx][:,1],fieldline_flf.fieldlines[idx][:,2],line_width=0.01,color=(0/255,0/255,0/255),tube_radius=0.002)


if fieldline_flf.geometry==1:
	for idx in range(len(fieldline_flf.fieldlines)):
		fieldline_X=fieldline_flf.fieldlines[idx][:,0]*np.sin(fieldline_flf.fieldlines[idx][:,1])*np.cos(fieldline_flf.fieldlines[idx][:,2])
		fieldline_Y=fieldline_flf.fieldlines[idx][:,0]*np.sin(fieldline_flf.fieldlines[idx][:,1])*np.sin(fieldline_flf.fieldlines[idx][:,2])
		fieldline_Z=fieldline_flf.fieldlines[idx][:,0]*np.cos(fieldline_flf.fieldlines[idx][:,1])
		mlab.plot3d(fieldline_X,fieldline_Y,fieldline_Z,line_width=0.01,color=(0/255,0/255,0/255),tube_radius=0.002)

mlab.show()

