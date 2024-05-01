

output_directory="./"

frames_per_step=15
frames_per_sec=60
pad_start_frames=0
pad_end_frames=0


import sys
UFiT_home_directory='/change/this/path'
sys.path.append(UFiT_home_directory)
from UFiT_Functions_Python import *
import matplotlib.pyplot as plt


Q_files = [f for f in os.listdir(output_directory) if os.path.isfile(os.path.join(output_directory, f)) and f.endswith('.flf')]
Q_files.sort()


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

fig=plt.figure("Emission",figsize=(16,16))

call_result=call(["mkdir","./anim_temp"])

max_Intensity=100.0*np.log(5000.0)/np.log(10.0)

plot_idx=0
for idx_f in range(len(Q_files)):

	fig.clf()

	Qfunc_full=np.zeros((len_vX,len_vY,len_vZ))

	UFiT_flf=read_UFiT_output(UFiTfile_name=os.path.join(output_directory,Q_files[idx_f]))

	for idx in range(len_indices):
		vr=np.sqrt(vX[indices[idx,0]]**2+vY[indices[idx,1]]**2+vZ[indices[idx,2]]**2)
		Qfunc_full[indices[idx,0],indices[idx,1],indices[idx,2]]=np.log(np.clip(abs(UFiT_flf.Q_out[idx]),1.0,5000.0))/np.log(10.0)/np.exp(vr)

	result=np.sum(Qfunc_full,axis=0)
	if idx_f==0:
		max_Intensity=max(result.flatten())

	ax1=fig.gca()
	plt.pcolormesh(vY,vZ,result.T,cmap="Greys_r",rasterized=True,vmin=0.0,vmax=max_Intensity)
	theta_c=np.linspace(0,2.0*np.pi,num=100)
	ax1.fill(1.0*np.cos(theta_c),1.0*np.sin(theta_c),color="#000000")
	ax1.set_xticks([])
	ax1.set_yticks([])

	plt.savefig("./anim_temp/img{:03d}.png".format(plot_idx), format="png", dpi=50,bbox_inches='tight',pad_inches=0.1)

	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(1,frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==len(Q_files)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1

call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","Emission_Q3D_Animate.mp4"])
call_result=call(["rm","-r","./anim_temp/"])

