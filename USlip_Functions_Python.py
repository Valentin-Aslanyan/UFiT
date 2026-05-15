

import numpy as np
import struct
import os
import ctypes


class USlip_call_input:	#input to call USlip from Python
	def __init__(self):
		self.str_mx=300
		self.geometry=0
		self.Bfile_type=-1
		self.grid_regular=True
		self.grid_separate=False
		self.periodic_X=False
		self.periodic_Y=False
		self.periodic_Z=False
		self.periodic_PHI=False
		self.print_devices=False
		self.num_proc=1
		self.B_filename="ufit.bin"
		self.out_directory=""
		self.load_B=False
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
		self.sz11=2
		self.sz12=2
		self.sz13=2
		self.sz21=2
		self.sz22=2
		self.sz23=2
		self.sz31=2
		self.sz32=2
		self.sz33=2
		self.grid1_1=np.zeros((2),dtype="double")
		self.grid1_2=np.zeros((2),dtype="double")
		self.grid1_3=np.zeros((2),dtype="double")
		self.grid2_1=np.zeros((2),dtype="double")
		self.grid2_2=np.zeros((2),dtype="double")
		self.grid2_3=np.zeros((2),dtype="double")
		self.grid3_1=np.zeros((2),dtype="double")
		self.grid3_2=np.zeros((2),dtype="double")
		self.grid3_3=np.zeros((2),dtype="double")
		self.B_grid1=np.zeros((2,2,2),dtype="double")
		self.B_grid2=np.zeros((2,2,2),dtype="double")
		self.B_grid3=np.zeros((2,2,2),dtype="double")
		self.write_output = False
		self.return_output = False
		self.j_grid = np.zeros((2,2,2,3),dtype="double")
		self.sigma_grid = np.zeros((2,2,2,3),dtype="double")
		self.sigmac_grid = np.zeros((2,2,2,3),dtype="double")
		self.j_grid_ir = np.zeros((1,2,2,2,3),dtype="double")
		self.sigma_grid_ir = np.zeros((1,2,2,2,3),dtype="double")
		self.sigmaalpha_grid_ir = np.zeros((1,2,2,2,3),dtype="double")


def call_USlip(shared_lib_path,USlip_input):
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
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.c_int, 
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
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.c_int, 
				ctypes.c_int, 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double)
				]

	if USlip_input.grid_regular==True or USlip_input.grid_regular==1:
		if USlip_input.grid_separate==True or USlip_input.grid_separate==1:
			USlip_input.sz11=len(USlip_input.grid1_1)
			USlip_input.sz12=len(USlip_input.grid1_2)
			USlip_input.sz13=len(USlip_input.grid1_3)
			USlip_input.sz21=len(USlip_input.grid2_1)
			USlip_input.sz22=len(USlip_input.grid2_2)
			USlip_input.sz23=len(USlip_input.grid2_3)
			USlip_input.sz31=len(USlip_input.grid3_1)
			USlip_input.sz32=len(USlip_input.grid3_2)
			USlip_input.sz33=len(USlip_input.grid3_3)
		else:
			USlip_input.sz1=len(USlip_input.grid1)
			USlip_input.sz2=len(USlip_input.grid2)
			USlip_input.sz3=len(USlip_input.grid3)
	else:
		if USlip_input.grid_separate==True or USlip_input.grid_separate==1:
			print("Irregular and separate grids incompatible - for now")
			return USlip_input
		else:
			USlip_input.num_blocks=np.shape(USlip_input.B_grid_ir)[4]
			USlip_input.sz1=np.shape(USlip_input.B_grid_ir)[1]
			USlip_input.sz2=np.shape(USlip_input.B_grid_ir)[2]
			USlip_input.sz3=np.shape(USlip_input.B_grid_ir)[3]

	USlip_input.grid1=np.array(USlip_input.grid1,dtype="double")
	USlip_input.grid2=np.array(USlip_input.grid2,dtype="double")
	USlip_input.grid3=np.array(USlip_input.grid3,dtype="double")
	USlip_input.B_grid=np.array(USlip_input.B_grid,dtype="double")
	USlip_input.grid1_ir=np.array(USlip_input.grid1_ir,dtype="double")
	USlip_input.grid2_ir=np.array(USlip_input.grid2_ir,dtype="double")
	USlip_input.grid3_ir=np.array(USlip_input.grid3_ir,dtype="double")
	USlip_input.B_grid_ir=np.array(USlip_input.B_grid_ir,dtype="double")
	USlip_input.grid1_1=np.array(USlip_input.grid1_1,dtype="double")
	USlip_input.grid1_2=np.array(USlip_input.grid1_2,dtype="double")
	USlip_input.grid1_3=np.array(USlip_input.grid1_3,dtype="double")
	USlip_input.grid2_1=np.array(USlip_input.grid2_1,dtype="double")
	USlip_input.grid2_2=np.array(USlip_input.grid2_2,dtype="double")
	USlip_input.grid2_3=np.array(USlip_input.grid2_3,dtype="double")
	USlip_input.grid3_1=np.array(USlip_input.grid3_1,dtype="double")
	USlip_input.grid3_2=np.array(USlip_input.grid3_2,dtype="double")
	USlip_input.grid3_3=np.array(USlip_input.grid3_3,dtype="double")
	USlip_input.B_grid1=np.array(USlip_input.B_grid1,dtype="double")
	USlip_input.B_grid2=np.array(USlip_input.B_grid2,dtype="double")
	USlip_input.B_grid3=np.array(USlip_input.B_grid3,dtype="double")
	if USlip_input.grid_regular==True or USlip_input.grid_regular==1:
		USlip_input.j_grid = np.zeros((USlip_input.sz3,USlip_input.sz2,USlip_input.sz1,3),dtype="double")
		USlip_input.sigma_grid = np.zeros((USlip_input.sz3,USlip_input.sz2,USlip_input.sz1,3),dtype="double")
		USlip_input.sigmaalpha_grid = np.zeros((USlip_input.sz3,USlip_input.sz2,USlip_input.sz1,3),dtype="double")
	else:
		USlip_input.j_grid_ir = np.zeros((USlip_input.num_blocks,USlip_input.sz3,USlip_input.sz2,USlip_input.sz1,3),dtype="double")
		USlip_input.sigma_grid_ir = np.zeros((USlip_input.num_blocks,USlip_input.sz3,USlip_input.sz2,USlip_input.sz1,3),dtype="double")
		USlip_input.sigmaalpha_grid_ir = np.zeros((USlip_input.num_blocks,USlip_input.sz3,USlip_input.sz2,USlip_input.sz1,3),dtype="double")


	fortranlib.USlip_Python_Callable(ctypes.c_int(USlip_input.geometry), 
					ctypes.c_int(USlip_input.Bfile_type), 
					ctypes.c_int(bool_to_int(USlip_input.grid_regular)), 
					ctypes.c_int(bool_to_int(USlip_input.grid_separate)), 
					ctypes.c_int(bool_to_int(USlip_input.periodic_X)), 
					ctypes.c_int(bool_to_int(USlip_input.periodic_Y)), 
					ctypes.c_int(bool_to_int(USlip_input.periodic_Z)), 
					ctypes.c_int(bool_to_int(USlip_input.periodic_PHI)), 
					ctypes.c_int(bool_to_int(USlip_input.print_devices)), 
					ctypes.c_int(USlip_input.num_proc), 
					convert_string(USlip_input.B_filename,USlip_input.str_mx), 
					convert_string(USlip_input.out_directory,USlip_input.str_mx), 
					ctypes.c_int(bool_to_int(USlip_input.load_B)), 
					ctypes.c_int(USlip_input.sz1), 
					ctypes.c_int(USlip_input.sz2), 
					ctypes.c_int(USlip_input.sz3), 
					ctypes.c_int(USlip_input.num_blocks), 
					USlip_input.grid1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid3.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(USlip_input.B_grid.transpose(3,2,1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(USlip_input.grid1_ir.transpose(1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(USlip_input.grid2_ir.transpose(1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(USlip_input.grid3_ir.transpose(1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(USlip_input.B_grid_ir.transpose(4,3,2,1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					ctypes.c_int(USlip_input.sz11), 
					ctypes.c_int(USlip_input.sz12), 
					ctypes.c_int(USlip_input.sz13), 
					ctypes.c_int(USlip_input.sz21), 
					ctypes.c_int(USlip_input.sz22), 
					ctypes.c_int(USlip_input.sz23), 
					ctypes.c_int(USlip_input.sz31), 
					ctypes.c_int(USlip_input.sz32), 
					ctypes.c_int(USlip_input.sz33), 
					USlip_input.grid1_1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid1_2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid1_3.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid2_1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid2_2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid2_3.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid3_1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid3_2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.grid3_3.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(USlip_input.B_grid1.transpose(2,1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(USlip_input.B_grid2.transpose(2,1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					np.ascontiguousarray(USlip_input.B_grid3.transpose(2,1,0)).ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					ctypes.c_int(bool_to_int(USlip_input.write_output)), 
					ctypes.c_int(bool_to_int(USlip_input.return_output)),
					USlip_input.j_grid.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.sigma_grid.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.sigmaalpha_grid.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.j_grid_ir.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.sigma_grid_ir.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					USlip_input.sigmaalpha_grid_ir.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					)

	USlip_input.j_grid=USlip_input.j_grid.transpose(3,2,1,0)
	USlip_input.sigma_grid=USlip_input.sigma_grid.transpose(3,2,1,0)
	USlip_input.sigmaalpha_grid=USlip_input.sigmaalpha_grid.transpose(3,2,1,0)
	USlip_input.j_grid_ir=USlip_input.j_grid_ir.transpose(4,3,2,1,0)
	USlip_input.sigma_grid_ir=USlip_input.sigma_grid_ir.transpose(4,3,2,1,0)
	USlip_input.sigmaalpha_grid_ir=USlip_input.sigmaalpha_grid_ir.transpose(4,3,2,1,0)

	return USlip_input
	
	



