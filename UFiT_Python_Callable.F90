!To compile without the makefile, run the following line (without leading exclamation mark):
!gfortran -shared -fPIC UFiT_Functions_Fortran.F90 Run_UFiT.F90 -O3 -fopenmp -o UFiT_Python_Callable.so

subroutine UFiT_Python_Callable(IN_geometry, &
                                IN_Bfile_type, &
                                IN_input_type, &
                                IN_grid_regular, &
                                IN_grid_separate, &
                                IN_periodic_X, &
                                IN_periodic_Y, &
                                IN_periodic_Z, &
                                IN_periodic_PHI, &
                                IN_print_devices, &
                                IN_save_endpoints, &
                                IN_save_Q, &
                                IN_save_fieldlines, &
                                IN_save_connection, &
                                IN_user_defined, &
                                IN_check_starts, &
                                IN_normalized_B, &
                                IN_include_curvature, &
                                IN_num_proc, &
                                IN_MAX_STEPS, &
                                IN_step_size, &
                                IN_B_filename, &
                                IN_out_filename, &
                                IN_load_B, &
                                IN_numin_tot, &
                                IN_numin1, &
                                IN_numin2, &
                                IN_numin3, &
                                IN_coord1, &
                                IN_coord2, &
                                IN_coord3, &
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
                                OUT_fieldline_endpoints, &
                                OUT_fieldline_allpos, &
                                OUT_fieldline_pts, &
                                OUT_fieldline_ptn, &
                                OUT_fieldline_Q, &
                                OUT_fieldline_connection, &
                                OUT_fieldline_user)  bind(c, name='UFiT_Python_Callable')

      use iso_c_binding
      use UFiT_Functions_Fortran
      implicit none

      integer(c_int), intent(in), value :: IN_geometry
      integer(c_int), intent(in), value :: IN_Bfile_type
      integer(c_int), intent(in), value :: IN_input_type
      integer(c_int), intent(in), value :: IN_grid_regular
      integer(c_int), intent(in), value :: IN_grid_separate
      integer(c_int), intent(in), value :: IN_periodic_X
      integer(c_int), intent(in), value :: IN_periodic_Y
      integer(c_int), intent(in), value :: IN_periodic_Z
      integer(c_int), intent(in), value :: IN_periodic_PHI
      integer(c_int), intent(in), value :: IN_print_devices
      integer(c_int), intent(in), value :: IN_save_endpoints
      integer(c_int), intent(in), value :: IN_save_Q
      integer(c_int), intent(in), value :: IN_save_fieldlines
      integer(c_int), intent(in), value :: IN_save_connection
      integer(c_int), intent(in), value :: IN_user_defined
      integer(c_int), intent(in), value :: IN_check_starts
      integer(c_int), intent(in), value :: IN_normalized_B
      integer(c_int), intent(in), value :: IN_include_curvature
      integer(c_int), intent(in), value :: IN_num_proc
      integer(c_int), intent(in), value :: IN_MAX_STEPS
      REAL(c_double), intent(in), value :: IN_step_size
      CHARACTER(kind=C_CHAR), intent(in) ::  IN_B_filename(str_mx)
      CHARACTER(kind=C_CHAR), intent(in) ::  IN_out_filename(str_mx)
      integer(c_int), intent(in), value :: IN_load_B
      integer(c_int), intent(in), value :: IN_numin_tot
      integer(c_int), intent(in), value :: IN_numin1
      integer(c_int), intent(in), value :: IN_numin2
      integer(c_int), intent(in), value :: IN_numin3
      REAL(c_double), intent(in) :: IN_coord1(IN_numin1)
      REAL(c_double), intent(in) :: IN_coord2(IN_numin2)
      REAL(c_double), intent(in) :: IN_coord3(IN_numin3)
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
      REAL(c_double), intent(out) :: OUT_fieldline_endpoints(6,IN_numin_tot)
      REAL(c_double), intent(out) :: OUT_fieldline_allpos(3,2*IN_MAX_STEPS+1,IN_numin_tot)
      integer(c_int), intent(out) :: OUT_fieldline_pts(IN_numin_tot)
      integer(c_int), intent(out) :: OUT_fieldline_ptn(IN_numin_tot)
      REAL(c_double), intent(out) :: OUT_fieldline_Q(IN_numin_tot)
      BYTE, intent(out) :: OUT_fieldline_connection(IN_numin_tot)
      REAL(c_double), intent(out) :: OUT_fieldline_user(1,IN_numin_tot)


      geometry = IN_geometry
      Bfile_type = IN_Bfile_type
      input_type = IN_input_type
      grid_regular = cint_to_logical(IN_grid_regular)
      grid_separate = cint_to_logical(IN_grid_separate)
      periodic_X = cint_to_logical(IN_periodic_X)
      periodic_Y = cint_to_logical(IN_periodic_Y)
      periodic_Z = cint_to_logical(IN_periodic_Z)
      periodic_PHI = cint_to_logical(IN_periodic_PHI)
      print_devices = cint_to_logical(IN_print_devices)
      save_endpoints = cint_to_logical(IN_save_endpoints)
      save_Q = cint_to_logical(IN_save_Q)
      save_fieldlines = cint_to_logical(IN_save_fieldlines)
      save_connection = cint_to_logical(IN_save_connection)
      user_defined = cint_to_logical(IN_user_defined)
      check_starts = cint_to_logical(IN_check_starts)
      normalized_B = cint_to_logical(IN_normalized_B)
      include_curvature = cint_to_logical(IN_include_curvature)
      num_proc = IN_num_proc
      MAX_STEPS = IN_MAX_STEPS
      step_size = IN_step_size
      B_filename = TRIM(transfer(IN_B_filename(1:str_mx), B_filename))
      out_filename = TRIM(transfer(IN_out_filename(1:str_mx), out_filename))

      numin1 = IN_numin1
      numin2 = IN_numin2
      numin3 = IN_numin3
      numin_tot = IN_numin_tot
      if (input_type .eq. 0) then
        ALLOCATE(coord1_in(numin_tot))
        ALLOCATE(coord2_in(numin_tot))
        ALLOCATE(coord3_in(numin_tot))
        coord1_in(:) = IN_coord1(:)
        coord2_in(:) = IN_coord2(:)
        coord3_in(:) = IN_coord3(:)
      else
        ALLOCATE(coord1_in(numin1))
        ALLOCATE(coord2_in(numin2))
        ALLOCATE(coord3_in(numin3))
        coord1_in(:) = IN_coord1(:)
        coord2_in(:) = IN_coord2(:)
        coord3_in(:) = IN_coord3(:)
      end if

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
      call get_available_resource
      call run_trace
      if (IN_write_output .eq. 1) then
        call write_output
      end if
      if (IN_return_output .eq. 1) then
        if (save_endpoints) then
          OUT_fieldline_endpoints(:,:) = fieldline_endpoints(:,:)
        end if
        if (save_fieldlines) then
          OUT_fieldline_allpos(:,:,:) = fieldline_allpos(:,:,:)
          OUT_fieldline_pts(:) = fieldline_pts(:)
          OUT_fieldline_ptn(:) = fieldline_ptn(:)
        end if
        if (save_Q) then
          OUT_fieldline_Q(:) = fieldline_Q(:)
        end if
        if (save_connection) then
          OUT_fieldline_connection(:) = fieldline_connection(:)
        end if
        if (user_defined) then
          OUT_fieldline_user(:,:) = fieldline_user(:,:)
        end if
      end if

      call cleanup

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

end subroutine UFiT_Python_Callable

