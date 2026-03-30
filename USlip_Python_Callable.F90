
subroutine USlip_Python_Callable(IN_geometry, &
                                 IN_Bfile_type, &
                                 IN_grid_regular, &
                                 IN_grid_separate, &
                                 IN_periodic_X, &
                                 IN_periodic_Y, &
                                 IN_periodic_Z, &
                                 IN_periodic_PHI, &
                                 IN_print_devices, &
                                 IN_num_proc, &
                                 IN_B_filename, &
                                 IN_out_directory, &
                                 IN_load_B, &
                                 IN_sz1, &
                                 IN_sz2, &
                                 IN_sz3, &
                                 IN_num_blocks, &
                                 IN_grid1, &
                                 IN_grid2, &
                                 IN_grid3, &
                                 IN_B_grid, &
                                 IN_grid1_ir, &
                                 IN_grid2_ir, &
                                 IN_grid3_ir, &
                                 IN_B_grid_ir, &
                                 IN_sz11, &
                                 IN_sz12, &
                                 IN_sz13, &
                                 IN_sz21, &
                                 IN_sz22, &
                                 IN_sz23, &
                                 IN_sz31, &
                                 IN_sz32, &
                                 IN_sz33, &
                                 IN_grid1_1, &
                                 IN_grid1_2, &
                                 IN_grid1_3, &
                                 IN_grid2_1, &
                                 IN_grid2_2, &
                                 IN_grid2_3, &
                                 IN_grid3_1, &
                                 IN_grid3_2, &
                                 IN_grid3_3, &
                                 IN_B_grid1, &
                                 IN_B_grid2, &
                                 IN_B_grid3, &
                                 IN_write_output, &
                                 IN_return_output, &
                                 OUT_j_grid, &
                                 OUT_sigma_grid, &
                                 OUT_sigmafac_grid, &
                                 OUT_j_grid_ir, &
                                 OUT_sigma_grid_ir, &
                                 OUT_sigmafac_grid_ir)  bind(c, name='USlip_Python_Callable')


      use iso_c_binding
      use UFiT_Functions_Fortran
      use USlip_Functions_Fortran
      implicit none

      integer(c_int), intent(in), value :: IN_geometry
      integer(c_int), intent(in), value :: IN_Bfile_type
      integer(c_int), intent(in), value :: IN_grid_regular
      integer(c_int), intent(in), value :: IN_grid_separate
      integer(c_int), intent(in), value :: IN_periodic_X
      integer(c_int), intent(in), value :: IN_periodic_Y
      integer(c_int), intent(in), value :: IN_periodic_Z
      integer(c_int), intent(in), value :: IN_periodic_PHI
      integer(c_int), intent(in), value :: IN_print_devices
      integer(c_int), intent(in), value :: IN_num_proc
      CHARACTER(kind=C_CHAR), intent(in) ::  IN_B_filename(str_mx)
      CHARACTER(kind=C_CHAR), intent(in) ::  IN_out_directory(str_mx)
      integer(c_int), intent(in), value :: IN_load_B
      integer(c_int), intent(in), value :: IN_sz1
      integer(c_int), intent(in), value :: IN_sz2
      integer(c_int), intent(in), value :: IN_sz3
      integer(c_int), intent(in), value :: IN_num_blocks
      REAL(c_double), intent(in) :: IN_grid1(IN_sz1)
      REAL(c_double), intent(in) :: IN_grid2(IN_sz2)
      REAL(c_double), intent(in) :: IN_grid3(IN_sz3)
      REAL(c_double), intent(in) :: IN_B_grid(3,IN_sz1,IN_sz2,IN_sz3)
      REAL(c_double), intent(in) :: IN_grid1_ir(2,IN_num_blocks)
      REAL(c_double), intent(in) :: IN_grid2_ir(2,IN_num_blocks)
      REAL(c_double), intent(in) :: IN_grid3_ir(2,IN_num_blocks)
      REAL(c_double), intent(in) :: IN_B_grid_ir(3,IN_sz1,IN_sz2,IN_sz3,IN_num_blocks)
      integer(c_int), intent(in), value :: IN_sz11
      integer(c_int), intent(in), value :: IN_sz21
      integer(c_int), intent(in), value :: IN_sz31
      integer(c_int), intent(in), value :: IN_sz12
      integer(c_int), intent(in), value :: IN_sz22
      integer(c_int), intent(in), value :: IN_sz32
      integer(c_int), intent(in), value :: IN_sz13
      integer(c_int), intent(in), value :: IN_sz23
      integer(c_int), intent(in), value :: IN_sz33
      REAL(c_double), intent(in) :: IN_grid1_1(IN_sz11)
      REAL(c_double), intent(in) :: IN_grid1_2(IN_sz12)
      REAL(c_double), intent(in) :: IN_grid1_3(IN_sz13)
      REAL(c_double), intent(in) :: IN_grid2_1(IN_sz21)
      REAL(c_double), intent(in) :: IN_grid2_2(IN_sz22)
      REAL(c_double), intent(in) :: IN_grid2_3(IN_sz23)
      REAL(c_double), intent(in) :: IN_grid3_1(IN_sz31)
      REAL(c_double), intent(in) :: IN_grid3_2(IN_sz32)
      REAL(c_double), intent(in) :: IN_grid3_3(IN_sz33)
      REAL(c_double), intent(in) :: IN_B_grid1(IN_sz11,IN_sz12,IN_sz13)
      REAL(c_double), intent(in) :: IN_B_grid2(IN_sz21,IN_sz22,IN_sz23)
      REAL(c_double), intent(in) :: IN_B_grid3(IN_sz31,IN_sz32,IN_sz33)
      integer(c_int), intent(in), value :: IN_write_output
      integer(c_int), intent(in), value :: IN_return_output
      REAL(c_double), intent(out) :: OUT_j_grid(3,IN_sz1,IN_sz2,IN_sz3)
      REAL(c_double), intent(out) :: OUT_sigma_grid(3,IN_sz1,IN_sz2,IN_sz3)
      REAL(c_double), intent(out) :: OUT_sigmafac_grid(3,IN_sz1,IN_sz2,IN_sz3)
      REAL(c_double), intent(out) :: OUT_j_grid_ir(3,IN_sz1,IN_sz2,IN_sz3,IN_num_blocks)
      REAL(c_double), intent(out) :: OUT_sigma_grid_ir(3,IN_sz1,IN_sz2,IN_sz3,IN_num_blocks)
      REAL(c_double), intent(out) :: OUT_sigmafac_grid_ir(3,IN_sz1,IN_sz2,IN_sz3,IN_num_blocks)


      geometry = IN_geometry
      Bfile_type = IN_Bfile_type
      grid_regular = cint_to_logical(IN_grid_regular)
      grid_separate = cint_to_logical(IN_grid_separate)
      periodic_X = cint_to_logical(IN_periodic_X)
      periodic_Y = cint_to_logical(IN_periodic_Y)
      periodic_Z = cint_to_logical(IN_periodic_Z)
      periodic_PHI = cint_to_logical(IN_periodic_PHI)
      print_devices = cint_to_logical(IN_print_devices)
      num_proc = IN_num_proc
      B_filename = TRIM(transfer(IN_B_filename(1:str_mx), B_filename))
      out_directory = TRIM(transfer(IN_out_directory(1:str_mx), out_directory))

      if (IN_load_B .eq. 1) then
        call load_Bfield
      else
        if (grid_separate) then
          sz_11 = IN_sz11
          sz_12 = IN_sz12
          sz_13 = IN_sz13
          sz_21 = IN_sz21
          sz_22 = IN_sz22
          sz_23 = IN_sz23
          sz_31 = IN_sz31
          sz_32 = IN_sz32
          sz_33 = IN_sz33
          if (grid_regular) then
            ALLOCATE(grid1_1(sz_11))
            ALLOCATE(grid1_2(sz_12))
            ALLOCATE(grid1_3(sz_13))
            ALLOCATE(grid2_1(sz_21))
            ALLOCATE(grid2_2(sz_22))
            ALLOCATE(grid2_3(sz_23))
            ALLOCATE(grid3_1(sz_31))
            ALLOCATE(grid3_2(sz_32))
            ALLOCATE(grid3_3(sz_33))
            ALLOCATE(B_grid1(sz_11,sz_12,sz_13))
            ALLOCATE(B_grid2(sz_21,sz_22,sz_23))
            ALLOCATE(B_grid3(sz_31,sz_32,sz_33))
            grid1_1(:) = IN_grid1_1(:)
            grid1_2(:) = IN_grid1_2(:)
            grid1_3(:) = IN_grid1_3(:)
            grid2_1(:) = IN_grid2_1(:)
            grid2_2(:) = IN_grid2_2(:)
            grid2_3(:) = IN_grid2_3(:)
            grid3_1(:) = IN_grid3_1(:)
            grid3_2(:) = IN_grid3_2(:)
            grid3_3(:) = IN_grid3_3(:)
            B_grid1(:,:,:) = IN_B_grid1(:,:,:)
            B_grid2(:,:,:) = IN_B_grid2(:,:,:)
            B_grid3(:,:,:) = IN_B_grid3(:,:,:)
          else
            print *, 'Irregular and separate grids incompatible - for now'
          end if
        else !grid not separate
          sz_1 = IN_sz1
          sz_2 = IN_sz2
          sz_3 = IN_sz3
          if (grid_regular) then
            ALLOCATE(grid1(sz_1))
            ALLOCATE(grid2(sz_2))
            ALLOCATE(grid3(sz_3))
            ALLOCATE(B_grid(3,sz_1,sz_2,sz_3))
            grid1(:) = IN_grid1(:)
            grid2(:) = IN_grid2(:)
            grid3(:) = IN_grid3(:)
            B_grid(:,:,:,:) = IN_B_grid(:,:,:,:)
          else
            num_blocks = IN_num_blocks
            ALLOCATE(grid1_ir(2,num_blocks))
            ALLOCATE(grid2_ir(2,num_blocks))
            ALLOCATE(grid3_ir(2,num_blocks))
            ALLOCATE(B_grid_ir(3,sz_1,sz_2,sz_3,num_blocks))
            grid1_ir(:,:) = IN_grid1_ir(:,:)
            grid2_ir(:,:) = IN_grid1_ir(:,:)
            grid3_ir(:,:) = IN_grid1_ir(:,:)
            B_grid_ir(:,:,:,:,:) = IN_B_grid_ir(:,:,:,:,:)
          end if
        end if
      end if
      call process_Bfield
      call process_Bfield_slip
      call get_available_resource
      call slip_calculation
      if (IN_write_output .eq. 1) then
        call write_slip_output
      end if
      if (IN_return_output .eq. 1) then
        if (grid_regular) then
          OUT_j_grid(:,:,:,:) = j_grid(:,:,:,:)
          OUT_sigma_grid(:,:,:,:) = sigma_grid(:,:,:,:)
          OUT_sigmafac_grid(:,:,:,:) = sigmafac_grid(:,:,:,:)
        else
          OUT_j_grid_ir(:,:,:,:,:) = j_grid_ir(:,:,:,:,:)
          OUT_sigma_grid_ir(:,:,:,:,:) = sigma_grid_ir(:,:,:,:,:)
          OUT_sigmafac_grid_ir(:,:,:,:,:) = sigmafac_grid_ir(:,:,:,:,:)
        end if
      end if

      call cleanup
      call slip_cleanup

      contains



      function cint_to_logical(int_in)

        integer(c_int) :: int_in

        LOGICAL :: cint_to_logical

        if (int_in .eq. 1) then
          cint_to_logical = .true.
        else
          cint_to_logical = .false.
        end if

      end function cint_to_logical

end subroutine USlip_Python_Callable

