
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
      call cleanup

end program Run_UFiT

