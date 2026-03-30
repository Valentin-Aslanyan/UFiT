
program Run_USlip

      use UFiT_Functions_Fortran
      use USlip_Functions_Fortran
      implicit none

      call initialize_variables
      call parse_slip_command_args
      call parse_slip_command_file
      call load_Bfield
      call process_Bfield
      call process_Bfield_slip
      call get_available_resource
      call slip_calculation
      call write_slip_output
      call cleanup
      call slip_cleanup

end program Run_USlip

