




import sys
UFiT_home_directory='/change/this/path'
sys.path.append(UFiT_home_directory)
from UFiT_Functions_Python import *
import matplotlib.pyplot as plt


idxfile=open('ufit_emission.idx','rb')
len_vX=struct.unpack('i',idxfile.read(4))[0]
vX=np.array(struct.unpack('d'*len_vX,idxfile.read(8*len_vX)))
len_vY=struct.unpack('i',idxfile.read(4))[0]
vY=np.array(struct.unpack('d'*len_vY,idxfile.read(8*len_vY)))
len_vZ=struct.unpack('i',idxfile.read(4))[0]
vZ=np.array(struct.unpack('d'*len_vZ,idxfile.read(8*len_vZ)))
len_indices=struct.unpack('i',idxfile.read(4))[0]
indices=np.zeros((len_indices,3),dtype=np.int32)
for idx in range(len_indices):
	indices[idx,:]=struct.unpack('iii',idxfile.read(12))[:]
idxfile.close()

Qfunc_full=np.zeros((len_vX,len_vY,len_vZ))

UFiT_flf=read_UFiT_output(UFiTfile_name='ufit_emission000.flf')

for idx in range(len_indices):
	vr=np.sqrt(vX[indices[idx,0]]**2+vY[indices[idx,1]]**2+vZ[indices[idx,2]]**2)
	Qfunc_full[indices[idx,0],indices[idx,1],indices[idx,2]]=np.log(np.clip(abs(UFiT_flf.Q_out[idx]),1.0,5000.0))/np.log(10.0)/np.exp(vr)

result=np.sum(Qfunc_full,axis=0)

fig=plt.figure("Emission",figsize=(8,8))
ax1=fig.gca()
plt.pcolormesh(vY,vZ,result.T,cmap="Greys_r",rasterized=True)
theta_c=np.linspace(0,2.0*np.pi,num=100)
ax1.fill(1.0*np.cos(theta_c),1.0*np.sin(theta_c),color="#000000")
ax1.set_xticks([])
ax1.set_yticks([])

plt.savefig("Emission_Q3D.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)
plt.show()
