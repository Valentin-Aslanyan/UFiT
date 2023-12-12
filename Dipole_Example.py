

"""
Calculate Q in the 1D, axisymmetric dipole case
Compare ordinary and normalized B mode
"""


R_min=1.01
R_ss=2.5
num_r=61
num_t=1810
num_p=4

numin=1048
in_th_start=0.85
in_th_end=0.9

plot_fieldlines=False


import sys
UFiT_home_directory='/change/this/path'
sys.path.append(UFiT_home_directory)
from UFiT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def Q_analytic_axisymmetric(R,theta,R0,Rss):
	C=0.5*np.log((2.0*Rss**3+R**3)/(3.0*R*Rss**2))
	theta_SL=np.arcsin(np.exp(-C))
	e2C=np.exp(2.0*C)
	sin2th=(np.sin(theta))**2
	cos2th=(np.cos(theta))**2
	Q=(Rss/R)**2*np.exp(2.0*C)*(1.0+(cos2th/(1.0-e2C*sin2th))*(1.0+((125.0/8.0/R**3-1.0)**2*sin2th)/((125.0/4.0/R**3+1.0)**2*cos2th)))*((125.0/4.0/Rss**3+1.0)*np.sqrt(1.0-e2C*sin2th))/np.sqrt((125.0/4.0/R**3+1.0)**2*cos2th+(125.0/8.0/R**3-1.0)**2*sin2th)
	Q[(theta>theta_SL) & (theta<np.pi-theta_SL)]=2.0
	return Q
	


grid_r=np.linspace(R_min,R_ss,num=num_r)
grid_t=np.linspace(0.0,np.pi,num=num_t)
grid_p=np.linspace(0.0,2.0*np.pi,num=num_p)

B=np.zeros((3,num_r,num_t,num_p))
B[0,:,:,:]=((R_min**3/grid_r**3)*(2.0*R_ss**3+grid_r**3)/(2.0*R_ss**3+R_min**3))[:,None,None]*np.cos(grid_t)[None,:,None]*np.ones((num_p))[None,None,:]
B[1,:,:,:]=((R_min**3/grid_r**3)*(R_ss**3-grid_r**3)/(2.0*R_ss**3+R_min**3))[:,None,None]*np.sin(grid_t)[None,:,None]*np.ones((num_p))[None,None,:]


r_in=[R_min]
theta_in=np.linspace(in_th_start,in_th_end,num=numin)
phi_in=[0.5] #doesn't matter
write_input_regular(r_in,theta_in,phi_in,infile_name='ufitdipole.inp')

Q_analytic=Q_analytic_axisymmetric(R_min,theta_in,R_min,R_ss)


write_B_file(grid_r,grid_t,grid_p,B,Bfile_name='ufit_dipole.bin')
os.system("./Run_UFiT -dl 0.005 -ms 5000 -np 4 -pp -se -sf -sq -g 1 -b ufit_dipole.bin -i ufitdipole.inp -o ufit_dipole.flf")
os.system("./Run_UFiT -dl 0.005 -ms 5000 -np 4 -pp -nb -se -sf -sq -g 1 -b ufit_dipole.bin -i ufitdipole.inp -o ufit_dipole_norm.flf")
os.system("./Run_UFiT -dl 0.005 -ms 5000 -np 4 -pp -ic -se -sf -sq -g 1 -b ufit_dipole.bin -i ufitdipole.inp -o ufit_dipole_ic.flf")
os.system("./Run_UFiT -dl 0.005 -ms 5000 -np 4 -pp -ic -nb -se -sf -sq -g 1 -b ufit_dipole.bin -i ufitdipole.inp -o ufit_dipole_ic_norm.flf")

UFiT_flf=read_UFiT_output(UFiTfile_name='ufit_dipole.flf')
UFiT_flf_norm=read_UFiT_output(UFiTfile_name='ufit_dipole_norm.flf')
UFiT_flf_ic=read_UFiT_output(UFiTfile_name='ufit_dipole_ic.flf')
UFiT_flf_ic_norm=read_UFiT_output(UFiTfile_name='ufit_dipole_ic_norm.flf')

plt.figure("Q, no curvature term",figsize=(10,7))
plt.plot(UFiT_flf.coord2,np.log(np.clip(abs(UFiT_flf.Q_out),2,1E10))/np.log(10.0),color="black",label="$B$")
plt.plot(UFiT_flf_norm.coord2,np.log(np.clip(abs(UFiT_flf_norm.Q_out),2,1E10))/np.log(10.0),color="red",label="$\\hat{B}$")

plt.ylabel(r"$\log_{10}(Q)$",fontsize=20)
plt.xlabel(r"$\theta$ [$^{c}$]",fontsize=20)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.legend(frameon=False,fontsize=20)

plt.figure("Q, with curvature term",figsize=(10,7))
plt.plot(UFiT_flf_ic.coord2,np.log(np.clip(abs(UFiT_flf_ic.Q_out),2,1E10))/np.log(10.0),color="blue",label="IC $B$")
plt.plot(UFiT_flf_ic_norm.coord2,np.log(np.clip(abs(UFiT_flf_ic_norm.Q_out),2,1E10))/np.log(10.0),color="orange",label="IC $\\hat{B}$")
plt.plot(theta_in,np.log(np.clip(abs(Q_analytic),2,1E10))/np.log(10.0),'--',color="red",label="Q analytic")

plt.ylabel(r"$\log_{10}(Q)$",fontsize=20)
plt.xlabel(r"$\theta$ [$^{c}$]",fontsize=20)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.legend(frameon=False,fontsize=20)

if plot_fieldlines:
	plt.figure("fieldlines",figsize=(4,8))
	for idx in range(len(UFiT_flf.fieldlines)):
		fieldline_X=UFiT_flf.fieldlines[idx][:,0]*np.sin(UFiT_flf.fieldlines[idx][:,1])
		fieldline_Z=UFiT_flf.fieldlines[idx][:,0]*np.cos(UFiT_flf.fieldlines[idx][:,1])
		plt.plot(fieldline_X,fieldline_Z,color="grey")
	plot_theta=np.linspace(0.0,np.pi,num=100)
	plt.plot(R_min*np.sin(plot_theta),R_min*np.cos(plot_theta),color="black")
	plt.plot(R_ss*np.sin(plot_theta),R_ss*np.cos(plot_theta),color="black")
	plt.ylabel(r"$Z$",fontsize=20)
	plt.xlabel(r"$X$",fontsize=20)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)

plt.show()









