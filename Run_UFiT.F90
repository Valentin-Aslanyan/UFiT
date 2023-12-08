!To compile without the makefile, run the following line (without leading exclamation mark):
!gfortran UFiT_Functions_Fortran.F90 Run_UFiT.F90 -O3 -fopenmp -o Run_UFiT

program Run_UFiT

      use UFiT_Functions_Fortran
      implicit none

      call initialize_variables
      call parse_command_args
      call parse_command_file
      call load_input
      call load_Bfield
      call process_Bfield
      call get_available_resource
      call run_trace
      call write_output

end program Run_UFiT

