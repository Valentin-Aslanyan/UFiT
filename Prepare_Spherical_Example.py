

"""
Create an analytic B field grid and fieldline starts for a UFiT run
"""


import sys
UFiT_home_directory='/change/this/path'
sys.path.append(UFiT_home_directory)
from UFiT_Functions_Python import *


R_min=1.0
R_ss=2.5
num_r=61
num_t=181
num_p=361

numin_r=1
numin_t=300
numin_p=600
r_in_limits=[1.0,1.0]
theta_in_limits=[0.0,np.pi]
phi_in_limits=[0.0,2.0*np.pi]


#Create B-field grid
grid_r=np.linspace(R_min,R_ss,num=num_r)
grid_t=np.linspace(0.0,np.pi,num=num_t)
grid_p=np.linspace(0.0,2.0*np.pi,num=num_p)

#PFSS dipole
B=np.zeros((3,num_r,num_t,num_p))
dip_m0=10.0
B[0,:,:,:]=dip_m0*((R_min**3/grid_r**3)*(2.0*R_ss**3+grid_r**3)/(2.0*R_ss**3+R_min**3))[:,None,None]*np.cos(grid_t)[None,:,None]*np.ones((num_p))[None,None,:]
B[1,:,:,:]=dip_m0*((R_min**3/grid_r**3)*(R_ss**3-grid_r**3)/(2.0*R_ss**3+R_min**3))[:,None,None]*np.sin(grid_t)[None,:,None]*np.ones((num_p))[None,None,:]

#Complicated expression for dipole fields: convert to Cartesians, then back
grid_r3=grid_r[:,None,None]*np.ones((num_t))[None,:,None]*np.ones((num_p))[None,None,:]
grid_t3=np.ones((num_r))[:,None,None]*grid_t[None,:,None]*np.ones((num_p))[None,None,:]
grid_p3=np.ones((num_r))[:,None,None]*np.ones((num_t))[None,:,None]*grid_p[None,None,:]
grid_st3=np.sin(grid_t3)
grid_ct3=np.cos(grid_t3)
grid_sp3=np.sin(grid_p3)
grid_cp3=np.cos(grid_p3)
grid_x3=grid_r3*np.sin(grid_t3)*np.cos(grid_p3)
grid_y3=grid_r3*np.sin(grid_t3)*np.sin(grid_p3)
grid_z3=grid_r3*np.cos(grid_t3)
dip_m1=0.75*np.array([-np.sin(1.0*np.pi),np.cos(1.0*np.pi),0.0])
dip_r1=0.75*np.array([np.cos(1.0*np.pi),np.sin(1.0*np.pi),0.0])
vec_r2=(grid_x3-dip_r1[0])**2+(grid_y3-dip_r1[1])**2+(grid_z3-dip_r1[2])**2
vec_r3=vec_r2**1.5
m_dot_r=dip_m1[0]*grid_x3+dip_m1[1]*grid_y3+dip_m1[2]*grid_z3
B_x=3.0*(grid_x3-dip_r1[0])*m_dot_r/vec_r2/vec_r3-dip_m1[0]/vec_r3
B_y=3.0*(grid_y3-dip_r1[1])*m_dot_r/vec_r2/vec_r3-dip_m1[1]/vec_r3
B_z=3.0*(grid_z3-dip_r1[2])*m_dot_r/vec_r2/vec_r3-dip_m1[2]/vec_r3
B[0,:,:,:]+=grid_st3*grid_cp3*B_x+grid_st3*grid_sp3*B_y+grid_ct3*B_z
B[1,:,:,:]+=grid_ct3*grid_cp3*B_x+grid_ct3*grid_sp3*B_y+grid_cp3*B_z
B[2,:,:,:]+=-grid_sp3*B_x-grid_st3*B_y

#Finally, write the analytic B grid
write_B_file(grid_r,grid_t,grid_p,B,Bfile_name='Example.bin')


#Create file of fieldline starts
r_in=np.linspace(r_in_limits[0],r_in_limits[1],num=numin_r)
theta_in=np.linspace(theta_in_limits[0],theta_in_limits[1],num=numin_t)
phi_in=np.linspace(phi_in_limits[0],phi_in_limits[1],num=numin_p)


#Write a formatted file for the benefit of the user; you can write an unformatted one easily by uncommenting the following:
#write_input_regular(r_in,theta_in,phi_in,infile_name='Example.inp')
outfile=open("Example.inp","w")
print("Input type: 2",file=outfile)
print(numin_r,file=outfile)
print(numin_t,file=outfile)
print(numin_p,file=outfile)
print(str(r_in_limits[0])+','+str(r_in_limits[1]),file=outfile)
print(str(theta_in_limits[0])+','+str(theta_in_limits[1]),file=outfile)
print(str(phi_in_limits[0])+','+str(phi_in_limits[1]),file=outfile)
outfile.close()

#Smaller number of fieldlines to plot
outfile=open("Example2.inp","w")
print("Input type: 2",file=outfile)
print(1,file=outfile)
print(10,file=outfile)
print(10,file=outfile)
print("1.0,1.0",file=outfile)
print(str(0.3*np.pi)+','+str(0.7*np.pi),file=outfile)
print(str(0.8*np.pi)+','+str(1.2*np.pi),file=outfile)
outfile.close()














