

"""
Create a 2D plot of Q from a successful UFiT run
"""


import sys
UFiT_home_directory='/change/this/path'
sys.path.append(UFiT_home_directory)
from UFiT_Functions_Python import *
import matplotlib.pyplot as plt


UFiT_flf=read_UFiT_output(UFiTfile_name='Example.flf')
dim1,dim2,Q=get_Q_2D(UFiT_flf)
Q[np.isinf(Q)]=1E6
Q[np.isinf(-Q)]=-1E6
slogQ=np.sign(Q)*np.log(np.clip(abs(Q),2,1E6))/np.log(10.0)

plt.figure(figsize=(10,6))
if UFiT_flf.geometry==0:
	plt.pcolormesh(dim1,dim2,slogQ,cmap="bwr",vmin=-6,vmax=6,rasterized=True)
if UFiT_flf.geometry==1:
	plt.pcolormesh(dim2*180.0/np.pi,90.0-dim1*180.0/np.pi,np.transpose(slogQ),cmap="bwr",vmin=-6,vmax=6,rasterized=True)
cbar=plt.colorbar()
cbar.solids.set_rasterized(True)
cbar.ax.tick_params(labelsize=20,direction='in', bottom=True, top=True, left=True, right=True)
cbar.set_label(label="$\mathrm{slog}(Q)$",fontsize=20)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
plt.tick_params(axis='both', which='major',labelsize=20,direction='in',bottom=True, top=True, left=True, right=True)
plt.savefig("Spherical_Example.png", dpi=100,bbox_inches='tight',pad_inches=0.1)
plt.show()


