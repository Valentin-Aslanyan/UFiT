

integration_step=0.02
view_phi=0		#degrees
pixels_dim=(500,500)
R_max=2.5
target_B_file="ufit.bin"


import sys
sys.path.append('/change/this/path')
from UFiT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


R_min=1.0
view_phi_rad=view_phi*np.pi/180.0

#Integrate along (viewpoint) X, this axis points away from you if at theta=0,phi=0; will be rotated to view_phi
vX_raw=np.linspace(-R_max,R_max,num=int(2.0*R_max/integration_step),dtype=np.float64)
vY_raw=np.linspace(-R_max,R_max,num=pixels_dim[0],dtype=np.float64)
vZ_raw=np.linspace(-R_max,R_max,num=pixels_dim[1],dtype=np.float64)


r_actual=[]
t_actual=[]
p_actual=[]
indices_actual=[]
for idx_y in range(len(vY_raw)):
	print(idx_y," / ",len(vY_raw))
	for idx_z in range(len(vZ_raw)):
		for idx_x in range(len(vX_raw)):
			vr=np.sqrt(vX_raw[idx_x]**2+vY_raw[idx_y]**2+vZ_raw[idx_z]**2)
			vr_partial=np.sqrt(vY_raw[idx_y]**2+vZ_raw[idx_z]**2)
			if vr<=R_max and vr>=R_min:	#can add the following to manually occult center:  and vr_partial>=R_min
				vt=np.arccos(vZ_raw[idx_z]/vr)
				vp=(np.arctan2(vY_raw[idx_y],vX_raw[idx_x])+view_phi_rad) % (2.0*np.pi)
				r_actual.append(vr)
				t_actual.append(vt)
				p_actual.append(vp)
				indices_actual.append([idx_x,idx_y,idx_z])

r_actual=np.array(r_actual,dtype=np.float64)
t_actual=np.array(t_actual,dtype=np.float64)
p_actual=np.array(p_actual,dtype=np.float64)
indices_actual=np.array(indices_actual,dtype=np.int32)

write_input_irregular(r_actual,t_actual,p_actual,infile_name='ufit_emission.inp')
outfile=open('ufit_emission.idx','wb')
outfile.write(struct.pack('i',len(vX_raw)))
outfile.write(vX_raw)
outfile.write(struct.pack('i',len(vY_raw)))
outfile.write(vY_raw)
outfile.write(struct.pack('i',len(vZ_raw)))
outfile.write(vZ_raw)
outfile.write(struct.pack('i',np.shape(indices_actual)[0]))
for idx in range(np.shape(indices_actual)[0]):
	outfile.write(indices_actual[idx,:])
outfile.close()


print("./Run_UFiT -dl 0.005 -ms 10000 -np 4 -pp -nb -sq -g 1 -b "+target_B_file+" -i ufit_emission.inp -o ufit_emission.flf")
os.system("./Run_UFiT -dl 0.005 -ms 10000 -np 4 -pp -nb -sq -g 1 -b "+target_B_file+" -i ufit_emission.inp -o ufit_emission.flf")



