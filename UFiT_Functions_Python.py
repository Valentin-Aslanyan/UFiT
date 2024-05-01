

import numpy as np
import struct
import os
import ctypes


UFiT_float_type=np.float64


class UFiT_flf_output:	#fieldline file, reading result of UFiT
	def __init__(self):
		self.fltsz=8
		self.fltname='d'
		self.version='1.0'
		self.geometry=0
		self.input_type=1
		self.save_endpoints = False
		self.save_Q = False
		self.save_fieldlines = False
		self.save_connection = False
		self.user_defined = False
		self.periodic_X = False
		self.periodic_Y = False
		self.periodic_Z = False
		self.periodic_PHI = False
		self.numin_tot = 2
		self.numin1 = 2
		self.numin2 = 1
		self.numin3 = 1
		self.coord1 = np.zeros((2))
		self.coord2 = np.zeros((1))
		self.coord3 = np.zeros((1))
		self.endpoints = np.zeros((2,6))
		self.Q_out = np.zeros((2))
		self.fieldlines = []
		self.connection = np.zeros((2))
		self.num_ud_variables = 1
		self.fieldline_user = np.zeros((1,1),dtype="double")

	def get_Q_2D(self,slice_dim=None,slice_index=None):
		if self.save_Q==False:
			print("Q not saved, use -sq command line option")
			return np.zeros((1)),np.zeros((1)),np.zeros((1,1))
		elif self.input_type==0:
			print("Input type 0 (each fieldline start point is uniquely identified) is unsuitable for 2D grid; use -it 2")
			return np.zeros((1)),np.zeros((1)),np.zeros((1,1))
		elif slice_dim==None:
			if self.numin1==1:
				return self.coord2, self.coord3, self.Q_out.reshape(self.numin3,self.numin2)
			elif self.numin2==1:
				return self.coord1, self.coord3, self.Q_out.reshape(self.numin3,self.numin1)
			elif self.numin3==1:
				return self.coord1, self.coord2, self.Q_out.reshape(self.numin2,self.numin1)
			else:
				print("Input grid (fieldline starts) is 3D, must specify slice_dim and slice_index")
				return np.zeros((1)),np.zeros((1)),np.zeros((1,1))
		elif slice_index==None:
			print("Specify slice index")
			return np.zeros((1)),np.zeros((1)),np.zeros((1,1))
		elif slice_dim==0 or slice_dim.lower()=='x' or slice_dim.lower()=='r':
			return self.coord2, self.coord3, self.Q_out.reshape(self.numin3,self.numin2,self.numin1)[:,:,slice_index]
		elif slice_dim==1 or slice_dim.lower()=='y' or slice_dim.lower()[0]=='t':
			return self.coord1, self.coord3, self.Q_out.reshape(self.numin3,self.numin2,self.numin1)[:,slice_index,:]
		elif slice_dim==2 or slice_dim.lower()=='z' or slice_dim.lower()[0]=='p':
			return self.coord1, self.coord2, self.Q_out.reshape(self.numin3,self.numin2,self.numin1)[slice_index,:,:]
		else:
			print("slice_dim or slice_index incorrect")
			return np.zeros((1)),np.zeros((1)),np.zeros((1,1))


	def get_Q_3D(self):
		if self.save_Q==False:
			print("Q not saved, use -sq command line option")
			return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1))
		elif self.input_type==0:
			print("Input type 0 (each fieldline start point is uniquely identified) is unsuitable for 2D grid; use -it 2")
			return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1))
		else:
			return self.coord1, self.coord2, self.coord3, self.Q_out.reshape(self.numin3,self.numin2,self.numin1)


	def get_connection_3D(self):
		if self.save_connection==False:
			print("connection not saved, use -sc command line option")
			return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1))
		elif self.input_type==0:
			print("Input type 0 (each fieldline start point is uniquely identified) is unsuitable for 2D grid; use -it 2")
			return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1))
		else:
			return self.coord1, self.coord2, self.coord3, self.connection.reshape(self.numin3,self.numin2,self.numin1)


class UFiT_call_input:	#input to call UFiT from Python
	def __init__(self):
		self.str_mx=300
		self.geometry=0
		self.Bfile_type=-1
		self.input_type=1
		self.grid_regular=True
		self.periodic_X=False
		self.periodic_Y=False
		self.periodic_Z=False
		self.periodic_PHI=False
		self.print_devices=False
		self.save_endpoints=True
		self.save_Q=False
		self.save_fieldlines=False
		self.save_connection=False
		self.user_defined=False
		self.check_starts=False
		self.normalized_B=True
		self.include_curvature=False
		self.num_proc=1
		self.MAX_STEPS=5000
		self.step_size=0.005
		self.B_filename="ufit.bin"
		self.out_filename="ufit.flf"
		self.load_B=False
		self.numin_tot=2
		self.numin1=2
		self.numin2=1
		self.numin3=1
		self.coord1=np.zeros((1),dtype="double")
		self.coord2=np.zeros((1),dtype="double")
		self.coord3=np.zeros((1),dtype="double")
		self.sz1=2
		self.sz2=2
		self.sz3=2
		self.num_blocks=1
		self.grid1=np.zeros((2),dtype="double")
		self.grid2=np.zeros((2),dtype="double")
		self.grid3=np.zeros((2),dtype="double")
		self.B_grid=np.zeros((3,2,2,2),dtype="double")
		self.grid1_ir=np.zeros((2,1),dtype="double")
		self.grid2_ir=np.zeros((2,1),dtype="double")
		self.grid3_ir=np.zeros((2,1),dtype="double")
		self.B_grid_ir=np.zeros((3,2,2,2,1),dtype="double")
		self.write_output = False
		self.return_output = False
		self.endpoints = np.zeros((6,2),dtype="double")
		self.fieldline_allpos = np.zeros((1,10001,3),dtype="double")
		self.fieldline_pts = np.zeros((1),dtype="int32")
		self.fieldline_ptn = np.zeros((1),dtype="int32")
		self.Q_out = np.zeros((2),dtype="double")
		self.fieldline_connection = np.zeros((2),dtype="b")
		self.num_ud_variables=1
		self.fieldline_user = np.zeros((1,1),dtype="double")


def call_UFiT(shared_lib_path,UFiT_input):
	def bool_to_int(bool_in):
		if bool_in==True or bool_in==1:
			return 1
		else:
			return 0

	def convert_string(str_in,str_mx):
		string_bin=str_in.encode('utf-8').ljust(str_mx)
		return string_bin

	fortranlib = ctypes.CDLL(shared_lib_path)
	fortranlib.argtypes = [	ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_double, 
				type(b'A'), 
				type(b'A'), 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_int), 
				ctypes.POINTER(ctypes.c_int), 
				ctypes.POINTER(ctypes.c_double),
				ctypes.POINTER(ctypes.c_byte),
				ctypes.POINTER(ctypes.c_double)
				]

	if UFiT_input.grid_regular==True or UFiT_input.grid_regular==1:
		UFiT_input.sz1=len(UFiT_input.grid1)
		UFiT_input.sz2=len(UFiT_input.grid2)
		UFiT_input.sz3=len(UFiT_input.grid3)
	else:
		UFiT_input.num_blocks=np.shape(UFiT_input.B_grid_ir)[4]
		UFiT_input.sz1=np.shape(UFiT_input.B_grid_ir)[1]
		UFiT_input.sz2=np.shape(UFiT_input.B_grid_ir)[2]
		UFiT_input.sz3=np.shape(UFiT_input.B_grid_ir)[3]

	UFiT_input.coord1=np.array(UFiT_input.coord1,dtype="double")
	UFiT_input.coord2=np.array(UFiT_input.coord2,dtype="double")
	UFiT_input.coord3=np.array(UFiT_input.coord3,dtype="double")
	UFiT_input.grid1=np.array(UFiT_input.grid1,dtype="double")
	UFiT_input.grid2=np.array(UFiT_input.grid2,dtype="double")
	UFiT_input.grid3=np.array(UFiT_input.grid3,dtype="double")
	UFiT_input.B_grid=np.array(UFiT_input.B_grid,dtype="double")
	UFiT_input.grid1_ir=np.array(UFiT_input.grid1_ir,dtype="double")
	UFiT_input.grid2_ir=np.array(UFiT_input.grid2_ir,dtype="double")
	UFiT_input.grid3_ir=np.array(UFiT_input.grid3_ir,dtype="double")
	UFiT_input.B_grid_ir=np.array(UFiT_input.B_grid_ir,dtype="double")
	if UFiT_input.input_type==0:
		UFiT_input.numin_tot=min(len(UFiT_input.coord1),len(UFiT_input.coord2),len(UFiT_input.coord3))
		UFiT_input.numin1=UFiT_input.numin_tot
		UFiT_input.numin2=UFiT_input.numin_tot
		UFiT_input.numin3=UFiT_input.numin_tot
	else:
		UFiT_input.numin1=len(UFiT_input.coord1)
		UFiT_input.numin2=len(UFiT_input.coord2)
		UFiT_input.numin3=len(UFiT_input.coord3)
		UFiT_input.numin_tot=UFiT_input.numin1*UFiT_input.numin2*UFiT_input.numin3
	UFiT_input.endpoints = np.zeros((6,UFiT_input.numin_tot),dtype="double")
	if UFiT_input.save_fieldlines:
		UFiT_input.fieldline_allpos = np.zeros((UFiT_input.numin_tot,2*UFiT_input.MAX_STEPS+1,3),dtype="double")
		UFiT_input.fieldline_pts = np.zeros((UFiT_input.numin_tot),dtype="int32")
		UFiT_input.fieldline_ptn = np.zeros((UFiT_input.numin_tot),dtype="int32")
	if UFiT_input.save_Q:
		UFiT_input.Q_out = np.zeros((UFiT_input.numin_tot),dtype="double")
	if UFiT_input.save_connection:
		UFiT_input.fieldline_connection = np.zeros((UFiT_input.numin_tot),dtype="b")
	if UFiT_input.save_connection:
		UFiT_input.fieldline_user = np.zeros((1,UFiT_input.numin_tot),dtype="double")


	fortranlib.UFiT_Python_Callable(ctypes.c_int(UFiT_input.geometry), 
					ctypes.c_int(UFiT_input.Bfile_type), 
					ctypes.c_int(UFiT_input.input_type), 
					ctypes.c_int(bool_to_int(UFiT_input.grid_regular)), 
					ctypes.c_int(bool_to_int(UFiT_input.periodic_X)), 
					ctypes.c_int(bool_to_int(UFiT_input.periodic_Y)), 
					ctypes.c_int(bool_to_int(UFiT_input.periodic_Z)), 
					ctypes.c_int(bool_to_int(UFiT_input.periodic_PHI)), 
					ctypes.c_int(bool_to_int(UFiT_input.print_devices)), 
					ctypes.c_int(bool_to_int(UFiT_input.save_endpoints)), 
					ctypes.c_int(bool_to_int(UFiT_input.save_Q)), 
					ctypes.c_int(bool_to_int(UFiT_input.save_fieldlines)), 
					ctypes.c_int(bool_to_int(UFiT_input.save_connection)), 
					ctypes.c_int(bool_to_int(UFiT_input.user_defined)), 
					ctypes.c_int(bool_to_int(UFiT_input.check_starts)), 
					ctypes.c_int(bool_to_int(UFiT_input.normalized_B)), 
					ctypes.c_int(bool_to_int(UFiT_input.include_curvature)), 
					ctypes.c_int(UFiT_input.num_proc), 
					ctypes.c_int(UFiT_input.MAX_STEPS), 
					ctypes.c_double(UFiT_input.step_size), 
					convert_string(UFiT_input.B_filename,UFiT_input.str_mx), 
					convert_string(UFiT_input.out_filename,UFiT_input.str_mx), 
					ctypes.c_int(bool_to_int(UFiT_input.load_B)), 
					ctypes.c_int(UFiT_input.numin_tot), 
					ctypes.c_int(UFiT_input.numin1), 
					ctypes.c_int(UFiT_input.numin2), 
					ctypes.c_int(UFiT_input.numin3), 
					UFiT_input.coord1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					UFiT_input.coord2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					UFiT_input.coord3.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					ctypes.c_int(UFiT_input.sz1), 
					ctypes.c_int(UFiT_input.sz2), 
					ctypes.c_int(UFiT_input.sz3), 
					ctypes.c_int(UFiT_input.num_blocks), 
					UFiT_input.grid1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					UFiT_input.grid2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					UFiT_input.grid3.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(UFiT_input.B_grid.transpose(3,2,1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(UFiT_input.grid1_ir.transpose(1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(UFiT_input.grid2_ir.transpose(1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(UFiT_input.grid3_ir.transpose(1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(UFiT_input.B_grid_ir.transpose(4,3,2,1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					ctypes.c_int(bool_to_int(UFiT_input.write_output)), 
					ctypes.c_int(bool_to_int(UFiT_input.return_output)),
					UFiT_input.endpoints.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					UFiT_input.fieldline_allpos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					UFiT_input.fieldline_pts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), 
					UFiT_input.fieldline_ptn.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), 
					UFiT_input.Q_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
					UFiT_input.fieldline_connection.ctypes.data_as(ctypes.POINTER(ctypes.c_byte)),
					UFiT_input.fieldline_user.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
					)

	return UFiT_input
	
	


def write_input_irregular(coord1,coord2,coord3,infile_name='ufit.inp'):
	coord1_arr=np.array(coord1)
	coord2_arr=np.array(coord2)
	coord3_arr=np.array(coord3)
	num_inp=min(len(coord1_arr),len(coord2_arr),len(coord3_arr))
	ufit_infile=open(infile_name,'wb')
	ufit_infile.write(struct.pack('i',0))
	ufit_infile.write(struct.pack('i',num_inp))
	for idx in range(num_inp):
		ufit_infile.write(coord1_arr[idx].astype(UFiT_float_type))
		ufit_infile.write(coord2_arr[idx].astype(UFiT_float_type))
		ufit_infile.write(coord3_arr[idx].astype(UFiT_float_type))
	ufit_infile.close()


def write_input_regular(coord1,coord2,coord3,infile_name='ufit.inp'):
	coord1_arr=np.array(coord1)
	coord2_arr=np.array(coord2)
	coord3_arr=np.array(coord3)
	ufit_infile=open(infile_name,'wb')
	ufit_infile.write(struct.pack('i',1))
	ufit_infile.write(struct.pack('i',len(coord1_arr)))
	ufit_infile.write(struct.pack('i',len(coord2_arr)))
	ufit_infile.write(struct.pack('i',len(coord3_arr)))
	ufit_infile.write(coord1_arr.astype(UFiT_float_type))
	ufit_infile.write(coord2_arr.astype(UFiT_float_type))
	ufit_infile.write(coord3_arr.astype(UFiT_float_type))
	ufit_infile.close()


#Default magnetic field file for UFiT:
#Ends in .bin
#INT(32) size of coordinate 1, n1 (X, r)
#INT(32) size of coordinate 2, n2 (Y, theta)
#INT(32) size of coordinate 3, n3 (Z, phi)
#FLOAT(64)*n1 ... coordinate 1 array
#FLOAT(64)*n1 ... coordinate 2 array
#FLOAT(64)*n1 ... coordinate 3 array
#FLOAT(64)*3*n1*n2*n3 ... B field ordered as follows: 
#first write each component, then increment dim1, then dim2, then dim3, e.g.:
#Bx(x=1,y=1,z=1),By(x=1,y=1,z=1),Bz(x=1,y=1,z=1),Bx(x=2,y=1,z=1),By(x=2,y=1,z=1) ...
def write_B_file(coord1,coord2,coord3,B_in,Bfile_name='ufit.bin'):
	if np.shape(B_in)!=(3,len(coord1),len(coord2),len(coord3)):
		print("Check B dimensions; should be (3,x,y,z)")
	else:
		ufit_Bfile=open(Bfile_name,'wb')
		ufit_Bfile.write(struct.pack('i',len(coord1)))
		ufit_Bfile.write(struct.pack('i',len(coord2)))
		ufit_Bfile.write(struct.pack('i',len(coord3)))
		ufit_Bfile.write(coord1.astype(UFiT_float_type))
		ufit_Bfile.write(coord2.astype(UFiT_float_type))
		ufit_Bfile.write(coord3.astype(UFiT_float_type))
		#reverse order for Fortran; set to UFiT_float_type
		ufit_Bfile.write(np.array(B_in,dtype=UFiT_float_type).transpose(3,2,1,0).flatten())
		ufit_Bfile.close()


def read_UFiT_output(UFiTfile_name='ufit.flf',print_header=False):
	f_o=UFiT_flf_output()

	ufit_file=open(UFiTfile_name,'rb')
	header_len_prelim=512
	header=ufit_file.read(header_len_prelim).decode("utf-8")
	header_len=int(header.split(';')[0].split(' ')[-1])
	header=header+ufit_file.read(header_len-header_len_prelim).decode("utf-8")
	if print_header:
		print(header)

	for idx in range(1,len(header.split(';'))):
		if len(header.split(';')[idx].split(' '))==3:
			command_code=header.split(';')[idx].split(' ')[1]
			command_value=header.split(';')[idx].split(' ')[2]
			if command_code.lower() == 'n:':
				f_o.fltsz=int(command_value)
			elif command_code.lower() == 'g:':
				f_o.geometry=int(command_value)
			elif command_code.lower() == 'it:':
				f_o.input_type=int(command_value)
			elif command_code.lower() == 'ic:':
				if command_value[0].lower()=='t':
					f_o.include_curvature=True
			elif command_code.lower() == 'se:':
				if command_value[0].lower()=='t':
					f_o.save_endpoints=True
			elif command_code.lower() == 'sq:':
				if command_value[0].lower()=='t':
					f_o.save_Q=True
			elif command_code.lower() == 'sf:':
				if command_value[0].lower()=='t':
					f_o.save_fieldlines=True
			elif command_code.lower() == 'sc:':
				if command_value[0].lower()=='t':
					f_o.save_connection=True
			elif command_code.lower() == 'ud:':
				if command_value[0].lower()=='t':
					f_o.user_defined=True
			elif command_code.lower() == 'px:':
				if command_value[0].lower()=='t':
					f_o.periodic_X=True
			elif command_code.lower() == 'py:':
				if command_value[0].lower()=='t':
					f_o.periodic_Y=True
			elif command_code.lower() == 'pz:':
				if command_value[0].lower()=='t':
					f_o.periodic_Z=True
			elif command_code.lower() == 'pp:':
				if command_value[0].lower()=='t':
					f_o.periodic_PHI=True

	if f_o.fltsz==8:
		f_o.fltname='d'
	else:
		f_o.fltname='f'

	if f_o.input_type==0:
		f_o.numin_tot=struct.unpack('i',ufit_file.read(4))[0]
		f_o.coord1=np.array(struct.unpack(f_o.fltname*f_o.numin_tot,ufit_file.read(f_o.fltsz*f_o.numin_tot)))
		f_o.coord2=np.array(struct.unpack(f_o.fltname*f_o.numin_tot,ufit_file.read(f_o.fltsz*f_o.numin_tot)))
		f_o.coord3=np.array(struct.unpack(f_o.fltname*f_o.numin_tot,ufit_file.read(f_o.fltsz*f_o.numin_tot)))
	else:
		f_o.numin1=struct.unpack('i',ufit_file.read(4))[0]
		f_o.numin2=struct.unpack('i',ufit_file.read(4))[0]
		f_o.numin3=struct.unpack('i',ufit_file.read(4))[0]
		f_o.numin_tot=f_o.numin1*f_o.numin2*f_o.numin3
		f_o.coord1=np.array(struct.unpack(f_o.fltname*f_o.numin1,ufit_file.read(f_o.fltsz*f_o.numin1)))
		f_o.coord2=np.array(struct.unpack(f_o.fltname*f_o.numin2,ufit_file.read(f_o.fltsz*f_o.numin2)))
		f_o.coord3=np.array(struct.unpack(f_o.fltname*f_o.numin3,ufit_file.read(f_o.fltsz*f_o.numin3)))
	
	if f_o.save_endpoints:
		f_o.endpoints=np.array(struct.unpack(f_o.fltname*6*f_o.numin_tot,ufit_file.read(f_o.fltsz*6*f_o.numin_tot))).reshape(f_o.numin_tot,6)
	if f_o.save_Q:
		f_o.Q_out=np.array(struct.unpack(f_o.fltname*f_o.numin_tot,ufit_file.read(f_o.fltsz*f_o.numin_tot)))
	if f_o.save_fieldlines:
		num_fl=np.array(struct.unpack('i'*f_o.numin_tot,ufit_file.read(4*f_o.numin_tot)),dtype=np.int32)
		for idx in range(f_o.numin_tot):
			f_o.fieldlines.append(np.array(struct.unpack(f_o.fltname*3*num_fl[idx],ufit_file.read(f_o.fltsz*3*num_fl[idx]))).reshape(num_fl[idx],3))
	if f_o.save_connection:
		f_o.connection=np.array(struct.unpack('b'*f_o.numin_tot,ufit_file.read(f_o.numin_tot)))
	if f_o.user_defined:
		f_o.num_ud_variables=struct.unpack('i',ufit_file.read(4))[0]
		f_o.fieldline_user=np.array(struct.unpack(f_o.fltname*f_o.num_ud_variables*f_o.numin_tot,ufit_file.read(f_o.fltsz*f_o.num_ud_variables*f_o.numin_tot))).reshape(f_o.numin_tot,f_o.num_ud_variables)
	ufit_file.close()
	return f_o


def get_Q_2D(flf_output,slice_dim=None,slice_index=None):
	if flf_output.save_Q==False:
		print("Q not saved, use -sq command line option")
		return np.zeros((1)),np.zeros((1)),np.zeros((1,1))
	elif flf_output.input_type==0:
		print("Input type 0 (each fieldline start point is uniquely identified) is unsuitable for 2D grid; use -it 2")
		return np.zeros((1)),np.zeros((1)),np.zeros((1,1))
	elif slice_dim==None:
		if flf_output.numin1==1:
			return flf_output.coord2, flf_output.coord3, flf_output.Q_out.reshape(flf_output.numin3,flf_output.numin2)
		elif flf_output.numin2==1:
			return flf_output.coord1, flf_output.coord3, flf_output.Q_out.reshape(flf_output.numin3,flf_output.numin1)
		elif flf_output.numin3==1:
			return flf_output.coord1, flf_output.coord2, flf_output.Q_out.reshape(flf_output.numin2,flf_output.numin1)
		else:
			print("Input grid (fieldline starts) is 3D, must specify slice_dim and slice_index")
			return np.zeros((1)),np.zeros((1)),np.zeros((1,1))
	elif slice_index==None:
		print("Specify slice index")
		return np.zeros((1)),np.zeros((1)),np.zeros((1,1))
	elif slice_dim==0 or slice_dim.lower()=='x' or slice_dim.lower()=='r':
		return flf_output.coord2, flf_output.coord3, flf_output.Q_out.reshape(flf_output.numin3,flf_output.numin2,flf_output.numin1)[:,:,slice_index]
	elif slice_dim==1 or slice_dim.lower()=='y' or slice_dim.lower()[0]=='t':
		return flf_output.coord1, flf_output.coord3, flf_output.Q_out.reshape(flf_output.numin3,flf_output.numin2,flf_output.numin1)[:,slice_index,:]
	elif slice_dim==2 or slice_dim.lower()=='z' or slice_dim.lower()[0]=='p':
		return flf_output.coord1, flf_output.coord2, flf_output.Q_out.reshape(flf_output.numin3,flf_output.numin2,flf_output.numin1)[slice_index,:,:]
	else:
		print("slice_dim or slice_index incorrect")
		return np.zeros((1)),np.zeros((1)),np.zeros((1,1))


def get_Q_3D(flf_output):
	if flf_output.save_Q==False:
		print("Q not saved, use -sq command line option")
		return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1))
	elif flf_output.input_type==0:
		print("Input type 0 (each fieldline start point is uniquely identified) is unsuitable for 3D grid; use -it 2")
		return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1))
	else:
		return flf_output.coord1, flf_output.coord2, flf_output.coord3, flf_output.Q_out.reshape(flf_output.numin3,flf_output.numin2,flf_output.numin1)


def get_connection_3D(flf_output):
	if flf_output.save_connection==False:
		print("connection not saved, use -sc command line option")
		return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1))
	elif flf_output.input_type==0:
		print("Input type 0 (each fieldline start point is uniquely identified) is unsuitable for 3D grid; use -it 2")
		return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1))
	else:
		return flf_output.coord1, flf_output.coord2, flf_output.coord3, flf_output.connection.reshape(flf_output.numin3,flf_output.numin2,flf_output.numin1)


def get_endpoints_3D(flf_output):
	if flf_output.save_endpoints==False:
		print("endpoints not saved, use -se command line option")
		return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1,6))
	elif flf_output.input_type==0:
		print("Input type 0 (each fieldline start point is uniquely identified) is unsuitable for 3D grid; use -it 2")
		return np.zeros((1)),np.zeros((1)),np.zeros((1)),np.zeros((1,1,1,6))
	else:
		return flf_output.coord1, flf_output.coord2, flf_output.coord3, flf_output.endpoints.reshape(flf_output.numin3,flf_output.numin2,flf_output.numin1,6)




