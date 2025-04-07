!!
!Coordinate indexing (note: Fortran starts with 1): 
!  1 = X for Cartesian, r for spherical
!  2 = Y for Cartesian, theta for spherical
!  3 = Z for Cartesian, phi for spherical
!B(3,grid1,grid2,grid3)
!For irregular (unstructured) grids, there is an addition num_blocks dimension

module UFiT_Functions_Fortran
      USE UFiT_Definitions_Fortran
      USE UFiT_User_Functions
      USE OMP_LIB
#ifdef USE_NC
      USE netcdf
#endif
      implicit none


      !public :: 


      contains


      subroutine initialize_variables

      !Geometries: 0 = Cartesian; 1 = spherical; 2 = cylindrical (TODO); 3 =? spherical exponential? (TODO)
        geometry = 0
      !Bfile file types: -1 = determine automatically; 
      !0 = UFiT default; 
      !10 = DUMFRIC spherical;
      !20 = Lare3d cartesian;
      !30 = ARMS flicks, geometry from header;
        Bfile_type = -1
        input_type = 0
        grid_regular = .true.
        grid_separate = .false.
        periodic_X = .false.
        periodic_Y = .false.
        periodic_Z = .false.
        periodic_PHI = .false.
        print_devices = .false.     !Print number of CPUs, GPUs available
        save_endpoints = .false.    !Positions of endpoints of fieldline (at boundary)
        save_Q = .false.            !Squashing factor
        save_fieldlines = .false.   !Full points along fieldline
        save_connection = .false.   !How the fieldline is connected; 0 = closed, 1 = open, 2 = disconnected
        user_defined = .false.      !TODO
        check_starts = .false.      !Loop over start points and check that they are inside the B_grid
        normalized_B = .false.      !Use B_hat for the calculation of Q
        include_curvature = .false. !Use full curvature formula for spherical coordinates
        num_proc = 1
        MAX_STEPS = 5000
        step_size = 0.005_num
        cmd_filename = 'ufit.dat'   !trace options, if not specified by command line
        B_filename   = 'ufit.bin'   !magnetic field
        in_filename  = 'ufit.inp'   !field line start points
        out_filename = 'ufit.flf'   !output with fieldline, endpoints, Q, ...

      end subroutine initialize_variables


      subroutine parse_command_args

        INTEGER :: idx, num_args, stat
        CHARACTER(len=str_mx) :: arg

        num_args = iargc()

        if (num_args .le. 0) then
          read_command_file = .true.
        else
          idx = 1
          DO WHILE (idx .le. num_args)
            CALL getarg(idx, arg)
            select case (arg)

              case ('-g', '--geometry')
                idx = idx + 1
                CALL getarg(idx, arg)
                read(arg,*,iostat=stat) geometry
                if (stat .ne. 0) then
                  print *, 'unrecognised geometry option: ', TRIM(arg)
                end if

              case ('-px', '--periodic_x', '--periodic_X')
                periodic_X = .true.

              case ('-py', '--periodic_y', '--periodic_Y')
                periodic_Y = .true.

              case ('-pz', '--periodic_z', '--periodic_Z')
                periodic_Z = .true.

              case ('-pp', '--periodic_phi', '--periodic_PHI')
                periodic_PHI = .true.

              case ('-gs', '--grid_staggered', '--grid_separate')
                grid_separate = .true.

              case ('-r', '--resources')
                print_devices = .true.

              case ('-se', '--save_endpoints')
                save_endpoints = .true.

              case ('-sq', '--save_Q')
                save_Q = .true.

              case ('-sf', '--save_fieldlines')
                save_fieldlines = .true.

              case ('-sc', '--save_connection')
                save_connection = .true.

              case ('-ud', '--user_defined')
                user_defined = .true.

              case ('-cs', '--check_starts')
                check_starts = .true.

              case ('-nb', '--normalize_B')
                normalized_B = .true.

              case ('-ic', '--include_curvature')
                !Note: this option is now deprecated and curvature is included by default
                include_curvature = .true.

              case ('-np', '--num_proc', '--NUM_PROC')
                idx = idx + 1
                CALL getarg(idx, arg)
                read(arg,*,iostat=stat) num_proc
                if (stat .ne. 0) then
                  print *, 'unrecognised step number: ', TRIM(arg)
                end if

              case ('-ms', '--max_steps', '--MAX_STEPS')
                idx = idx + 1
                CALL getarg(idx, arg)
                read(arg,*,iostat=stat) MAX_STEPS
                if (stat .ne. 0) then
                  print *, 'unrecognised step number: ', TRIM(arg)
                end if

              case ('-c', '--command_file')
                idx = idx + 1
                CALL getarg(idx, arg)
                cmd_filename = TRIM(arg)
                read_command_file = .true.

              case ('-i', '--input_file')
                idx = idx + 1
                CALL getarg(idx, arg)
                in_filename = TRIM(arg)

              case ('-b', '--B_file')
                idx = idx + 1
                CALL getarg(idx, arg)
                B_filename = TRIM(arg)

              case ('-bt', '--Bfile_type')
                idx = idx + 1
                CALL getarg(idx, arg)
                read(arg,*,iostat=stat) Bfile_type
                if (stat .ne. 0) then
                  print *, 'unrecognised B file type: ', TRIM(arg)
                end if

              case ('-dl', '--step_size')
                idx = idx + 1
                CALL getarg(idx, arg)
                read(arg,*,iostat=stat) step_size
                if (stat .ne. 0) then
                  print *, 'unrecognised step size: ', TRIM(arg)
                end if

              case ('-o', '--output_file')
                idx = idx + 1
                CALL getarg(idx, arg)
                out_filename = TRIM(arg)

              case default
                print *, 'unrecognised command-line option: ', TRIM(arg)

            end select
            idx = idx + 1
          END DO
        end if

      end subroutine parse_command_args


      subroutine parse_command_file

        LOGICAL :: cfile_exists, continue_read
        INTEGER :: stat, stat2, cmd_unit
        CHARACTER(len=str_mx) :: arg, arg2

        IF (read_command_file) THEN
          continue_read = .true.
          inquire(file=cmd_filename, exist=cfile_exists)
          IF (.not. cfile_exists) THEN
            print *, 'Unable to open command file'
            print *, 'Specify a valid command file, or pass command line arguments'
            print *, 'Attemped to read filename: '
            print *, TRIM(cmd_filename)
            call EXIT(110)
          END IF
          open(newunit=cmd_unit,file=cmd_filename,form="formatted")

          DO WHILE (continue_read)
            read(cmd_unit, *, iostat=stat) arg, arg2
            IF (stat .eq. 0) THEN
              select case (arg)

                case ('G:')
                  read(arg2,*,iostat=stat2) geometry
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised geometry option: ', TRIM(arg2)
                  end if

                case ('PX:')
                  read(arg2,*,iostat=stat2) periodic_X
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised periodic X option: ', TRIM(arg2)
                  end if

                case ('PY:')
                  read(arg2,*,iostat=stat2) periodic_Y
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised periodic Y option: ', TRIM(arg2)
                  end if

                case ('PZ:')
                  read(arg2,*,iostat=stat2) periodic_Z
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised periodic Z option: ', TRIM(arg2)
                  end if

                case ('PP:')
                  read(arg2,*,iostat=stat2) periodic_PHI
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised periodic PHI option: ', TRIM(arg2)
                  end if

                case ('GS:')
                  read(arg2,*,iostat=stat2) grid_separate
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised grid separate/staggered option: ', TRIM(arg2)
                  end if

                case ('R:')
                  read(arg2,*,iostat=stat2) print_devices
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised resource option: ', TRIM(arg2)
                  end if

                case ('SE:')
                  read(arg2,*,iostat=stat2) save_endpoints
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised save endpoints option: ', TRIM(arg2)
                  end if

                case ('SQ:')
                  read(arg2,*,iostat=stat2) save_Q
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised save Q option: ', TRIM(arg2)
                  end if

                case ('SF:')
                  read(arg2,*,iostat=stat2) save_fieldlines
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised save fieldlines option: ', TRIM(arg2)
                  end if

                case ('SC:')
                  read(arg2,*,iostat=stat2) save_connection
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised save connection option: ', TRIM(arg2)
                  end if

                case ('UD:')
                  read(arg2,*,iostat=stat2) user_defined
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised user defined option: ', TRIM(arg2)
                  end if

                case ('CS:')
                  read(arg2,*,iostat=stat2) check_starts
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised check starts option: ', TRIM(arg2)
                  end if

                case ('NB:')
                  read(arg2,*,iostat=stat2) normalized_B
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised B normalization option: ', TRIM(arg2)
                  end if

                case ('IC:')
                !Note: this option is now deprecated and curvature is included by default
                  read(arg2,*,iostat=stat2) include_curvature
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised curvature option: ', TRIM(arg2)
                  end if

                case ('NP:')
                  read(arg2,*,iostat=stat2) num_proc
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised number of processors: ', TRIM(arg2)
                  end if

                case ('MS:')
                  read(arg2,*,iostat=stat2) MAX_STEPS
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised maximum steps: ', TRIM(arg2)
                  end if

                case ('C:')
                  cmd_filename = TRIM(arg2)
                  read_command_file = .true.

                case ('I:')
                  in_filename = TRIM(arg2)

                case ('B:')
                  B_filename = TRIM(arg2)

                case ('BT:')
                  read(arg2,*,iostat=stat2) Bfile_type
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised B file type: ', TRIM(arg2)
                  end if

                case ('DL:')
                  read(arg2,*,iostat=stat2) step_size
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised step size: ', TRIM(arg2)
                  end if

                case ('O:')
                  out_filename = TRIM(arg2)

                case default
                  print *, 'unrecognised command file option: ', TRIM(arg)

              end select
            ELSE
              continue_read = .false.
            END IF
          END DO

          close(cmd_unit)
        END IF

      end subroutine parse_command_file


      subroutine get_available_resource
      !Find the number of CPUs/GPUs

        IF (print_devices) THEN
          print *, "Number of CPUs available: ", OMP_GET_MAX_THREADS()
        END IF

        IF (num_proc .gt. OMP_GET_MAX_THREADS()) THEN
          num_proc = OMP_GET_MAX_THREADS()
        ELSE IF (num_proc .lt. 1) THEN
          num_proc = 1
        END IF

      end subroutine get_available_resource


      subroutine load_input

        LOGICAL :: ifile_exists, continue_read
        INTEGER :: in_unit, stat, idx, idx2
        REAL(num) :: start_in_x, start_in_y, start_in_z, end_in_x, end_in_y, end_in_z
        CHARACTER(len=str_mx) :: arg

      !Read input file (field line seed coordinates)
      !Input types: 
      !   0 = each point is explicitly by set of coordinates (X,Y,Z), (R,Theta,Phi) etc
      !   1 = seed coordinates are a regular grid; each axis is specified explicitly
      !   2 = seed coordinates are an evenly-spaced grid, specified by the start and end

      !Unformatted (binary) format by default
        inquire(file=in_filename, exist=ifile_exists)
        IF (.not. ifile_exists) THEN
          print *, 'Unable to open input file'
          print *, 'This is required for the field line starts'
          print *, 'Attemped to read filename: '
          print *, TRIM(in_filename)
          call EXIT(111)
        END IF
        open(newunit=in_unit,file=in_filename,access='stream')
        read(in_unit) input_type
        if (input_type .eq. 0) then
          read(in_unit) numin_tot
          ALLOCATE(coord1_in(numin_tot))
          ALLOCATE(coord2_in(numin_tot))
          ALLOCATE(coord3_in(numin_tot))
          DO idx = 1, numin_tot
            read(in_unit) coord1_in(idx)
            read(in_unit) coord2_in(idx)
            read(in_unit) coord3_in(idx)
          END DO
          close(in_unit)
        else if ((input_type .eq. 1)) then
          read(in_unit) numin1
          read(in_unit) numin2
          read(in_unit) numin3
          numin_tot = numin1*numin2*numin3
          ALLOCATE(coord1_in(numin1))
          ALLOCATE(coord2_in(numin2))
          ALLOCATE(coord3_in(numin3))
          read(in_unit) coord1_in(:)
          read(in_unit) coord2_in(:)
          read(in_unit) coord3_in(:)
          close(in_unit)
        else if ((input_type .eq. 2)) then
          read(in_unit) numin1
          read(in_unit) numin2
          read(in_unit) numin3
          numin_tot = numin1*numin2*numin3
          ALLOCATE(coord1_in(numin1))
          ALLOCATE(coord2_in(numin2))
          ALLOCATE(coord3_in(numin3))
          read(in_unit) start_in_x
          read(in_unit) end_in_x
          read(in_unit) start_in_y
          read(in_unit) end_in_y
          read(in_unit) start_in_z
          read(in_unit) end_in_z
          IF (numin1 .eq. 1) THEN
            coord1_in(1) = start_in_x
          ELSE
            DO idx = 1, numin1
              coord1_in(idx) = start_in_x + (end_in_x-start_in_x)*(idx - 1)/(numin1 - 1)
            END DO
          END IF
          IF (numin2 .eq. 1) THEN
            coord2_in(1) = start_in_y
          ELSE
            DO idx = 1, numin2
              coord2_in(idx) = start_in_y + (end_in_y-start_in_y)*(idx - 1)/(numin2 - 1)
            END DO
          END IF
          IF (numin3 .eq. 1) THEN
            coord3_in(1) = start_in_z
          ELSE
            DO idx = 1, numin3
              coord3_in(idx) = start_in_z + (end_in_z-start_in_z)*(idx - 1)/(numin3 - 1)
            END DO
          END IF

      !Formatted text input file
        else 
          close(in_unit)
          open(newunit=in_unit,file=in_filename,form="formatted")
          read(in_unit, '(A)', iostat=stat) arg
          idx2 = 1
          continue_read = .true.
          DO idx = str_mx, 1, -1
            if (continue_read .and. (arg(idx:idx) .ne. ' ')) then
              continue_read = .false.
            else if ((.not. continue_read) .and. (arg(idx:idx) .eq. ' ')) then
              idx2 = idx
              exit
            end if
          END DO
          read(arg(idx2:), *) input_type
          if (input_type .eq. 0) then
            read(in_unit, *, iostat=stat) numin_tot
            ALLOCATE(coord1_in(numin_tot))
            ALLOCATE(coord2_in(numin_tot))
            ALLOCATE(coord3_in(numin_tot))
            DO idx = 1, numin_tot
              read(in_unit, *, iostat=stat) coord1_in(idx), coord2_in(idx), coord3_in(idx)
            END DO
          else if ((input_type .eq. 1)) then
            read(in_unit, *, iostat=stat) numin1
            read(in_unit, *, iostat=stat) numin2
            read(in_unit, *, iostat=stat) numin3
            numin_tot = numin1*numin2*numin3
            ALLOCATE(coord1_in(numin1))
            ALLOCATE(coord2_in(numin2))
            ALLOCATE(coord3_in(numin3))
            DO idx = 1, numin1
              read(in_unit, *, iostat=stat) coord1_in(idx)
            END DO
            DO idx = 1, numin2
              read(in_unit, *, iostat=stat) coord2_in(idx)
            END DO
            DO idx = 1, numin3
              read(in_unit, *, iostat=stat) coord3_in(idx)
            END DO
          else if ((input_type .eq. 2)) then
            read(in_unit, *, iostat=stat) numin1
            read(in_unit, *, iostat=stat) numin2
            read(in_unit, *, iostat=stat) numin3
            numin_tot = numin1*numin2*numin3
            ALLOCATE(coord1_in(numin1))
            ALLOCATE(coord2_in(numin2))
            ALLOCATE(coord3_in(numin3))
            read(in_unit, *, iostat=stat) start_in_x, end_in_x
            read(in_unit, *, iostat=stat) start_in_y, end_in_y
            read(in_unit, *, iostat=stat) start_in_z, end_in_z
            IF (numin1 .eq. 1) THEN
              coord1_in(1) = start_in_x
            ELSE
              DO idx = 1, numin1
                coord1_in(idx) = start_in_x + (end_in_x-start_in_x)*(idx - 1)/(numin1 - 1)
              END DO
            END IF
            IF (numin2 .eq. 1) THEN
              coord2_in(1) = start_in_y
            ELSE
              DO idx = 1, numin2
                coord2_in(idx) = start_in_y + (end_in_y-start_in_y)*(idx - 1)/(numin2 - 1)
              END DO
            END IF
            IF (numin3 .eq. 1) THEN
              coord3_in(1) = start_in_z
            ELSE
              DO idx = 1, numin3
                coord3_in(idx) = start_in_z + (end_in_z-start_in_z)*(idx - 1)/(numin3 - 1)
              END DO
            END IF
          end if
          close(in_unit)
        end if

      end subroutine load_input


      subroutine load_Bfield

        LOGICAL :: bfile_exists, stop_found
        INTEGER :: Bfile_type_actual, B_unit, stp_idx, file_size

      !Automatically detect file type based on extension
        stop_found = .false.
        if (Bfile_type .eq. -1) then
          stp_idx = str_mx
          DO WHILE ((stp_idx .gt. 0) .and. (.not. stop_found))
            if (B_filename(stp_idx:stp_idx) .eq. '.') then
              stop_found = .true.
            else
              stp_idx = stp_idx - 1
            end if
          END DO
          if (TRIM(B_filename(LEN(TRIM(B_filename))-3:)) .eq. '.bin') then
            Bfile_type_actual = 0
          else if (TRIM(B_filename(LEN(TRIM(B_filename))-2:)) .eq. '.nc') then
            Bfile_type_actual = 10
          else if (TRIM(B_filename(LEN(TRIM(B_filename))-3:)) .eq. '.sdf') then
            Bfile_type_actual = 20
          else if (stp_idx .ge. 7) then
            if (B_filename(stp_idx-6:stp_idx-1) .eq. 'flicks') then
              Bfile_type_actual = 30
            else if (B_filename(stp_idx-6:stp_idx-1) .eq. 'bfield') then
              Bfile_type_actual = 31
            else
              print *, 'unable to automatically identify Bfile type: ', Bfile_type
              call EXIT(102)
            end if
          else
            print *, 'unable to automatically identify Bfile type '
            print *, 'Attemped to read filename: '
            print *, TRIM(B_filename)
            call EXIT(102)
          end if
        else if ((Bfile_type .eq. 0) .or. (Bfile_type .eq. 10) .or. (Bfile_type .eq. 20) .or. &
                 (Bfile_type .eq. 30) .or. (Bfile_type .eq. 31)) then
          Bfile_type_actual = Bfile_type
        else
          print *, 'unrecognised Bfile type: ', Bfile_type
          call EXIT(101)
        end if

      !UFiT unformatted
        if (Bfile_type_actual .eq. 0) then
          grid_regular = .true.
          inquire(file=B_filename, exist=bfile_exists, SIZE=file_size)
          IF (.not. bfile_exists) THEN
            print *, 'Unable to open B field file'
            print *, 'Attemped to read filename: '
            print *, TRIM(B_filename)
            call EXIT(112)
          END IF
          IF (file_size .lt. 36) THEN
            print *, 'B field file is too small'
            print *, 'Attemped to read filename: '
            print *, TRIM(B_filename)
            print *, 'Detected size: ', file_size
            call EXIT(112)
          END IF
          open(newunit=B_unit,file=B_filename,access='stream')
          IF (grid_separate) THEN
            read(B_unit) sz_11
            read(B_unit) sz_12
            read(B_unit) sz_13
            read(B_unit) sz_21
            read(B_unit) sz_22
            read(B_unit) sz_23
            read(B_unit) sz_31
            read(B_unit) sz_32
            read(B_unit) sz_33
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
            read(B_unit) grid1_1
            read(B_unit) grid1_2
            read(B_unit) grid1_3
            read(B_unit) grid2_1
            read(B_unit) grid2_2
            read(B_unit) grid2_3
            read(B_unit) grid3_1
            read(B_unit) grid3_2
            read(B_unit) grid3_3
            read(B_unit) B_grid1
            read(B_unit) B_grid2
            read(B_unit) B_grid3
          ELSE
            read(B_unit) sz_1
            read(B_unit) sz_2
            read(B_unit) sz_3
            ALLOCATE(grid1(sz_1))
            ALLOCATE(grid2(sz_2))
            ALLOCATE(grid3(sz_3))
            ALLOCATE(B_grid(3,sz_1,sz_2,sz_3))
            read(B_unit) grid1
            read(B_unit) grid2
            read(B_unit) grid3
            read(B_unit) B_grid
          END IF
          close(B_unit)
      !DUMFRIC spherical (3D)
        else if (Bfile_type_actual .eq. 10) then
          grid_regular = .true.
          call load_DUMFRIC
      !Lare3d cartesian (3D)
        else if (Bfile_type_actual .eq. 20) then
          grid_regular = .true.
          IF (grid_separate) THEN
            print *, 'Separate/staggered grid setting incompatible with Lare3d output file'
            print *, 'Attempting to use ordinary grid instead'
            grid_separate = .false.
          END IF
          call load_Lare3d
      !ARMS flicks, geometry from header (3D)
        else if ((Bfile_type_actual .eq. 30) .or. (Bfile_type_actual .eq. 31)) then
          grid_regular = .false.
          IF (grid_separate) THEN
            print *, 'Separate/staggered grid setting not yet implemented for ARMS output file'
            print *, 'Attempting to use ordinary grid instead'
            grid_separate = .false.
          END IF
          IF (Bfile_type_actual .eq. 30) THEN
            call load_ARMS_flicks
          ELSE IF (Bfile_type_actual .eq. 31) THEN
            call load_ARMS_bfield
          END IF
        end if

      end subroutine load_Bfield


      subroutine process_Bfield

        INTEGER :: idx, idx_b, idx1, idx2, idx3
        REAL(num) :: coord_in1min, coord_in1max, coord_in2min, coord_in2max
        REAL(num) :: coord_in3min, coord_in3max
        REAL(num), DIMENSION(:), ALLOCATABLE :: grid_temp
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: B_temp

      !Reverse grid to make it increasing everywhere
        if (grid_separate) then
          if (grid_regular) then
          !Reverse B_grid1 and related coordinates
            if (grid1_1(2) .lt. grid1_1(1)) then
              ALLOCATE(grid_temp(sz_11))
              ALLOCATE(B_temp(sz_11,sz_12,sz_13))
              grid_temp(:) = grid1_1(:)
              B_temp(:,:,:) = B_grid1(:,:,:)
              do idx1 = 1, sz_11
                grid1_1(idx1) = grid_temp(sz_11+1-idx1)
                B_grid1(idx1,:,:) = B_temp(sz_11+1-idx1,:,:)
              end do
              DEALLOCATE(grid_temp)
              DEALLOCATE(B_temp)
            end if
            if (grid1_2(2) .lt. grid1_2(1)) then
              ALLOCATE(grid_temp(sz_12))
              ALLOCATE(B_temp(sz_11,sz_12,sz_13))
              grid_temp(:) = grid1_2(:)
              B_temp(:,:,:) = B_grid1(:,:,:)
              do idx2 = 1, sz_12
                grid1_2(idx2) = grid_temp(sz_12+1-idx2)
                B_grid1(:,idx2,:) = B_temp(:,sz_12+1-idx2,:)
              end do
              DEALLOCATE(grid_temp)
              DEALLOCATE(B_temp)
            end if
            if (grid1_3(2) .lt. grid1_3(1)) then
              ALLOCATE(grid_temp(sz_13))
              ALLOCATE(B_temp(sz_11,sz_12,sz_13))
              grid_temp(:) = grid1_3(:)
              B_temp(:,:,:) = B_grid1(:,:,:)
              do idx3 = 1, sz_13
                grid1_3(idx3) = grid_temp(sz_13+1-idx3)
                B_grid1(:,:,idx3) = B_temp(:,:,sz_13+1-idx3)
              end do
              DEALLOCATE(grid_temp)
              DEALLOCATE(B_temp)
            end if
          !Reverse B_grid2 and related coordinates
            if (grid2_1(2) .lt. grid2_1(1)) then
              ALLOCATE(grid_temp(sz_21))
              ALLOCATE(B_temp(sz_21,sz_22,sz_23))
              grid_temp(:) = grid2_1(:)
              B_temp(:,:,:) = B_grid2(:,:,:)
              do idx1 = 1, sz_21
                grid2_1(idx1) = grid_temp(sz_21+1-idx1)
                B_grid2(idx1,:,:) = B_temp(sz_21+1-idx1,:,:)
              end do
              DEALLOCATE(grid_temp)
              DEALLOCATE(B_temp)
            end if
            if (grid2_2(2) .lt. grid2_2(1)) then
              ALLOCATE(grid_temp(sz_22))
              ALLOCATE(B_temp(sz_21,sz_22,sz_23))
              grid_temp(:) = grid2_2(:)
              B_temp(:,:,:) = B_grid2(:,:,:)
              do idx2 = 1, sz_22
                grid2_2(idx2) = grid_temp(sz_22+1-idx2)
                B_grid2(:,idx2,:) = B_temp(:,sz_22+1-idx2,:)
              end do
              DEALLOCATE(grid_temp)
              DEALLOCATE(B_temp)
            end if
            if (grid2_3(2) .lt. grid2_3(1)) then
              ALLOCATE(grid_temp(sz_23))
              ALLOCATE(B_temp(sz_21,sz_22,sz_23))
              grid_temp(:) = grid2_3(:)
              B_temp(:,:,:) = B_grid2(:,:,:)
              do idx3 = 1, sz_23
                grid2_3(idx3) = grid_temp(sz_23+1-idx3)
                B_grid2(:,:,idx3) = B_temp(:,:,sz_23+1-idx3)
              end do
              DEALLOCATE(grid_temp)
              DEALLOCATE(B_temp)
            end if
          !Reverse B_grid3 and related coordinates
            if (grid3_1(2) .lt. grid3_1(1)) then
              ALLOCATE(grid_temp(sz_31))
              ALLOCATE(B_temp(sz_31,sz_32,sz_33))
              grid_temp(:) = grid3_1(:)
              B_temp(:,:,:) = B_grid3(:,:,:)
              do idx1 = 1, sz_31
                grid3_1(idx1) = grid_temp(sz_31+1-idx1)
                B_grid3(idx1,:,:) = B_temp(sz_31+1-idx1,:,:)
              end do
              DEALLOCATE(grid_temp)
              DEALLOCATE(B_temp)
            end if
            if (grid3_2(2) .lt. grid3_2(1)) then
              ALLOCATE(grid_temp(sz_32))
              ALLOCATE(B_temp(sz_31,sz_32,sz_33))
              grid_temp(:) = grid3_2(:)
              B_temp(:,:,:) = B_grid3(:,:,:)
              do idx2 = 1, sz_32
                grid3_2(idx2) = grid_temp(sz_32+1-idx2)
                B_grid3(:,idx2,:) = B_temp(:,sz_32+1-idx2,:)
              end do
              DEALLOCATE(grid_temp)
              DEALLOCATE(B_temp)
            end if
            if (grid2_3(2) .lt. grid2_3(1)) then
              ALLOCATE(grid_temp(sz_33))
              ALLOCATE(B_temp(sz_31,sz_32,sz_33))
              grid_temp(:) = grid3_3(:)
              B_temp(:,:,:) = B_grid3(:,:,:)
              do idx3 = 1, sz_33
                grid3_3(idx3) = grid_temp(sz_33+1-idx3)
                B_grid3(:,:,idx3) = B_temp(:,:,sz_33+1-idx3)
              end do
              DEALLOCATE(grid_temp)
              DEALLOCATE(B_temp)
            end if

            grid1min = MAX(grid1_1(1),grid2_1(1),grid3_1(1))
            grid1max = MIN(grid1_1(sz_11),grid2_1(sz_21),grid3_1(sz_31))
            grid2min = MAX(grid1_2(1),grid2_2(1),grid3_2(1))
            grid2max = MIN(grid1_2(sz_12),grid2_2(sz_22),grid3_2(sz_32))
            grid3min = MAX(grid1_3(1),grid2_3(1),grid3_3(1))
            grid3max = MIN(grid1_3(sz_13),grid2_3(sz_23),grid3_3(sz_33))
          else !Grid irregular
            !TODO
          end if
        else !grid not separate
          if (grid_regular) then
            if (grid1(2) .lt. grid1(1)) then
              ALLOCATE(grid_temp(sz_1))
              grid_temp(:) = grid1(:)
              do idx1 = 1, sz_1
                grid1(idx1) = grid_temp(sz_1+1-idx1)
              end do
              DEALLOCATE(grid_temp)
              ALLOCATE(B_temp(sz_1,sz_2,sz_3))
              do idx = 1, 3
                B_temp(:,:,:) = B_grid(idx,:,:,:)
                do idx1 = 1, sz_1
                  B_grid(idx,idx1,:,:) = B_temp(sz_1+1-idx1,:,:)
                end do
              end do
              DEALLOCATE(B_temp)
            end if
            if (grid2(2) .lt. grid2(1)) then
              ALLOCATE(grid_temp(sz_2))
              grid_temp(:) = grid2(:)
              do idx2 = 1, sz_2
                grid2(idx2) = grid_temp(sz_2+1-idx2)
              end do
              DEALLOCATE(grid_temp)
              ALLOCATE(B_temp(sz_1,sz_2,sz_3))
              do idx = 1, 3
                B_temp(:,:,:) = B_grid(idx,:,:,:)
                do idx2 = 1, sz_2
                  B_grid(idx,:,idx2,:) = B_temp(:,sz_2+1-idx2,:)
                end do
              end do
              DEALLOCATE(B_temp)
            end if
            if (grid3(2) .lt. grid3(1)) then
              ALLOCATE(grid_temp(sz_3))
              grid_temp(:) = grid3(:)
              do idx3 = 1, sz_3
                grid3(idx3) = grid_temp(sz_3+1-idx3)
              end do
              DEALLOCATE(grid_temp)
              ALLOCATE(B_temp(sz_1,sz_2,sz_3))
              do idx = 1, 3
                B_temp(:,:,:) = B_grid(idx,:,:,:)
                do idx3 = 1, sz_3
                  B_grid(idx,:,:,idx3) = B_temp(:,:,sz_3+1-idx3)
                end do
              end do
              DEALLOCATE(B_temp)
            end if

            grid1min = grid1(1)
            grid1max = grid1(sz_1)
            grid2min = grid2(1)
            grid2max = grid2(sz_2)
            grid3min = grid3(1)
            grid3max = grid3(sz_3)
          else !Grid irregular
            grid1min = MINVAL(grid1_ir(:,1))
            grid1max = MAXVAL(grid1_ir(:,1))
            grid2min = MINVAL(grid2_ir(:,1))
            grid2max = MAXVAL(grid2_ir(:,1))
            grid3min = MINVAL(grid3_ir(:,1))
            grid3max = MAXVAL(grid3_ir(:,1))
            do idx_b=1,num_blocks
              if (grid1_ir(2,idx_b) .lt. grid1_ir(1,idx_b)) then
                ALLOCATE(grid_temp(2))
                grid_temp(:) = grid1_ir(:,idx_b)
                grid1_ir(1,idx_b) = grid_temp(2)
                grid1_ir(2,idx_b) = grid_temp(1)
                DEALLOCATE(grid_temp)
                ALLOCATE(B_temp(sz_1,sz_2,sz_3))
                do idx = 1, 3
                  B_temp(:,:,:) = B_grid_ir(idx,:,:,:,idx_b)
                  do idx1 = 1, sz_1
                    B_grid_ir(idx,idx1,:,:,idx_b) = B_temp(sz_1+1-idx1,:,:)
                  end do
                end do
                DEALLOCATE(B_temp)
              end if
              if (grid2_ir(2,idx_b) .lt. grid2_ir(1,idx_b)) then
                ALLOCATE(grid_temp(2))
                grid_temp(:) = grid2_ir(:,idx_b)
                grid2_ir(1,idx_b) = grid_temp(2)
                grid2_ir(2,idx_b) = grid_temp(1)
                DEALLOCATE(grid_temp)
                ALLOCATE(B_temp(sz_1,sz_2,sz_3))
                do idx = 1, 3
                  B_temp(:,:,:) = B_grid_ir(idx,:,:,:,idx_b)
                  do idx2 = 1, sz_2
                    B_grid_ir(idx,:,idx2,:,idx_b) = B_temp(:,sz_2+1-idx2,:)
                  end do
                end do
                DEALLOCATE(B_temp)
              end if
              if (grid3_ir(2,idx_b) .lt. grid3_ir(1,idx_b)) then
                ALLOCATE(grid_temp(2))
                grid_temp(:) = grid3_ir(:,idx_b)
                grid3_ir(1,idx_b) = grid_temp(2)
                grid3_ir(2,idx_b) = grid_temp(1)
                DEALLOCATE(grid_temp)
                ALLOCATE(B_temp(sz_1,sz_2,sz_3))
                do idx = 1, 3
                  B_temp(:,:,:) = B_grid_ir(idx,:,:,:,idx_b)
                  do idx3 = 1, sz_3
                    B_grid_ir(idx,:,:,idx3,idx_b) = B_temp(:,:,sz_3+1-idx3)
                  end do
                end do
                DEALLOCATE(B_temp)
              end if

              grid1min = MIN(grid1min,grid1_ir(1,idx_b))
              grid1max = MAX(grid1max,grid1_ir(2,idx_b))
              grid2min = MIN(grid2min,grid2_ir(1,idx_b))
              grid2max = MAX(grid2max,grid2_ir(2,idx_b))
              grid3min = MIN(grid3min,grid3_ir(1,idx_b))
              grid3max = MAX(grid3max,grid3_ir(2,idx_b))
            end do
          end if
        end if

        coord_width(1) = grid1max-grid1min
        coord_width(2) = grid2max-grid2min
        coord_width(3) = grid3max-grid3min

        if (check_starts) then
          coord_in1min=MINVAL(coord1_in)
          coord_in1max=MAXVAL(coord1_in)
          coord_in2min=MINVAL(coord2_in)
          coord_in2max=MAXVAL(coord2_in)
          coord_in3min=MINVAL(coord3_in)
          coord_in3max=MAXVAL(coord3_in)
          print *, 'Input grid has the following extent:'
          if (geometry .eq. 0) then
            print *, 'Dimension 1 (X), min: ',grid1min, ' max: ',grid1max
            print *, 'Dimension 2 (Y), min: ',grid2min, ' max: ',grid2max
            print *, 'Dimension 3 (Z), min: ',grid3min, ' max: ',grid3max
            if ((coord_in1min .lt. grid1min) .or. (coord_in1max .gt. grid1max) .or. &
                (coord_in2min .lt. grid2min) .or. (coord_in2max .gt. grid2max) .or. &
                (coord_in3min .lt. grid3min) .or. (coord_in3max .gt. grid3max)) then
              print *, 'Error in start points, they have the following extent:'
              print *, 'Dimension 1 (X), min: ',coord_in1min, ' max: ',coord_in1max
              print *, 'Dimension 2 (Y), min: ',coord_in2min, ' max: ',coord_in2max
              print *, 'Dimension 3 (Z), min: ',coord_in3min, ' max: ',coord_in3max
              print *, 'Calculation will run, but expect some errors in outputs'
            else
               print *, 'Start points OK; all inside input grid'
            end if 
          else if (geometry .eq. 1) then
            print *, 'Dimension 1 (R), min: ',grid1min, ' max: ',grid1max
            print *, 'Dimension 2 (Theta), min: ',grid2min, ' max: ',grid2max
            print *, 'Dimension 3 (Phi), min: ',grid3min, ' max: ',grid3max
            if ((coord_in1min .lt. grid1min) .or. (coord_in1max .gt. grid1max) .or. &
                (coord_in2min .lt. grid2min) .or. (coord_in2max .gt. grid2max) .or. &
                (coord_in3min .lt. grid3min) .or. (coord_in3max .gt. grid3max)) then
              print *, 'Error in start points, they have the following extent:'
              print *, 'Dimension 1 (R), min: ',coord_in1min, ' max: ',coord_in1max
              print *, 'Dimension 2 (Theta), min: ',coord_in2min, ' max: ',coord_in2max
              print *, 'Dimension 3 (Phi), min: ',coord_in3min, ' max: ',coord_in3max
              print *, 'Calculation will run, but expect some errors in outputs'
            else
               print *, 'Start points OK; all inside input grid'
            end if 
          end if

        end if

      end subroutine process_Bfield


      subroutine load_DUMFRIC

#ifdef USE_NC
      !Compiling with NETCDF enabled

        LOGICAL :: bfile_exists
        INTEGER :: stat_nc, nc_id, var_id, ndims, sz_1_temp, sz_2_temp, sz_3_temp
        integer, dimension(nf90_max_var_dims) :: dimids
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: B_temp

      !Get dimensions, allocate and read grids
        inquire(file=B_filename, exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open B field file'
          print *, 'Attemped to read filename: '
          print *, TRIM(B_filename)
          call EXIT(115)
        END IF
        IF (grid_separate) THEN
          stat_nc = NF90_OPEN(TRIM(B_filename), NF90_NOWRITE, nc_id)
          stat_nc = NF90_INQ_VARID(nc_id, 'r', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_1_temp)
          sz_11 = sz_1_temp
          stat_nc = NF90_INQ_VARID(nc_id, 'rc', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_1_temp)
          sz_21 = sz_1_temp
          sz_31 = sz_1_temp
          stat_nc = NF90_INQ_VARID(nc_id, 'th', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_2_temp)
          sz_22 = sz_2_temp
          stat_nc = NF90_INQ_VARID(nc_id, 'thc', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_2_temp)
          sz_12 = sz_2_temp
          sz_32 = sz_2_temp
          stat_nc = NF90_INQ_VARID(nc_id, 'ph', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_3_temp)
          sz_33 = sz_3_temp
          stat_nc = NF90_INQ_VARID(nc_id, 'phc', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_3_temp)
          sz_13 = sz_3_temp
          sz_23 = sz_3_temp
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

      !Get grids and B field directly in each dimension
          stat_nc = NF90_INQ_VARID(nc_id, 'r', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, grid1_1)
          stat_nc = NF90_INQ_VARID(nc_id, 'rc', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, grid2_1)
          grid3_1(:) = grid2_1(:)
          stat_nc = NF90_INQ_VARID(nc_id, 'th', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, grid2_2)
          stat_nc = NF90_INQ_VARID(nc_id, 'thc', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, grid1_2)
          grid3_2(:) = grid1_2(:)
          stat_nc = NF90_INQ_VARID(nc_id, 'ph', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, grid3_3)
          stat_nc = NF90_INQ_VARID(nc_id, 'phc', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, grid1_3)
          grid2_3(:) = grid1_3(:)
          stat_nc = NF90_INQ_VARID(nc_id, 'br', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, B_grid1)
          stat_nc = NF90_INQ_VARID(nc_id, 'bth', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, B_grid2)
          stat_nc = NF90_INQ_VARID(nc_id, 'bph', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, B_grid3)

          stat_nc = NF90_CLOSE(nc_id)
        ELSE
          stat_nc = NF90_OPEN(TRIM(B_filename), NF90_NOWRITE, nc_id)
          stat_nc = NF90_INQ_VARID(nc_id, 'r', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_1)
          stat_nc = NF90_INQ_VARID(nc_id, 'th', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_2)
          stat_nc = NF90_INQ_VARID(nc_id, 'ph', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_3)
          ALLOCATE(grid1(sz_1))
          ALLOCATE(grid2(sz_2))
          ALLOCATE(grid3(sz_3))
          ALLOCATE(B_grid(3,sz_1,sz_2,sz_3))
          stat_nc = NF90_INQ_VARID(nc_id, 'r', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, grid1)
          stat_nc = NF90_INQ_VARID(nc_id, 'th', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, grid2)
          stat_nc = NF90_INQ_VARID(nc_id, 'ph', var_id)
          stat_nc = NF90_GET_VAR(nc_id, var_id, grid3)

      !Get B field and average 4 vertices in each respective dimension
          stat_nc = NF90_INQ_VARID(nc_id, 'br', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_1_temp)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(2), len = sz_2_temp)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(3), len = sz_3_temp)
          ALLOCATE(B_temp(sz_1_temp, sz_2_temp, sz_3_temp))
          stat_nc = NF90_GET_VAR(nc_id, var_id, B_temp)
          B_grid(1,:,:,:) = 0.25_num * (B_temp(:,2:,2:)+B_temp(:,2:,:sz_3_temp-1)+ &
                          B_temp(:,:sz_2_temp-1,2:)+B_temp(:,:sz_2_temp-1,:sz_3_temp-1))
          DEALLOCATE(B_temp)
          stat_nc = NF90_INQ_VARID(nc_id, 'bth', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_1_temp)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(2), len = sz_2_temp)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(3), len = sz_3_temp)
          ALLOCATE(B_temp(sz_1_temp, sz_2_temp, sz_3_temp))
          stat_nc = NF90_GET_VAR(nc_id, var_id, B_temp)
          B_grid(2,:,:,:) = 0.25_num * (B_temp(2:,:,2:)+B_temp(2:,:,:sz_3_temp-1)+ &
                          B_temp(:sz_1_temp-1,:,2:)+B_temp(:sz_1_temp-1,:,:sz_3_temp-1))
          DEALLOCATE(B_temp)
          stat_nc = NF90_INQ_VARID(nc_id, 'bph', var_id)
          stat_nc = nf90_inquire_variable(nc_id, var_id, ndims=ndims, dimids=dimids)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(1), len = sz_1_temp)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(2), len = sz_2_temp)
          stat_nc = nf90_inquire_dimension(nc_id, dimids(3), len = sz_3_temp)
          ALLOCATE(B_temp(sz_1_temp, sz_2_temp, sz_3_temp))
          stat_nc = NF90_GET_VAR(nc_id, var_id, B_temp)
          B_grid(3,:,:,:) = 0.25_num * (B_temp(2:,2:,:)+B_temp(2:,:sz_2_temp-1,:)+ &
                          B_temp(:sz_1_temp-1,2:,:)+B_temp(:sz_1_temp-1,:sz_2_temp-1,:))
          DEALLOCATE(B_temp)

          stat_nc = NF90_CLOSE(nc_id)
        END IF

#else
      !Default compilation: NETCDF disabled
        print *, 'Cannot read DUMFRIC file; UFiT was compiled without netCDF'
        print *, 'Ensure netCDF libraries are correct, then run: make USE_NC=True'
        call EXIT(121)
#endif

      end subroutine load_DUMFRIC


      subroutine load_Lare3d

        LOGICAL :: bfile_exists
        INTEGER :: idx_b
      !Lare3d/Warwick CFSA SDF variables
        INTEGER(1) :: restart_flag, subdomain_file
        INTEGER(4) :: L3D_unit, endianness, sdf_version, sdf_revision,summary_size, nblocks
        INTEGER(4) :: block_header_length, step, jobid1, jobid2, string_length, code_io_version
        INTEGER(4) :: blocktype, datatype, ndims, block_info_length, geometry_type
        INTEGER(8) :: first_block_location, next_block_location, current_block_location
        INTEGER(8) :: summary_location, data_location, data_length
        CHARACTER(len=4) :: sdf1_str
        CHARACTER(len=32) :: code_name
        CHARACTER(len=32) :: block_id
        CHARACTER(len=32) :: data_description
        REAL(8) :: time
        REAL(8), DIMENSION(:), ALLOCATABLE :: mults
        INTEGER(4), DIMENSION(:), ALLOCATABLE :: dims
        CHARACTER(len=:), ALLOCATABLE :: block_name
        REAL(4), DIMENSION(:), ALLOCATABLE :: arr_r4_1d
        REAL(8), DIMENSION(:), ALLOCATABLE :: arr_r8_1d
        REAL(16), DIMENSION(:), ALLOCATABLE :: arr_r16_1d
        INTEGER(4), DIMENSION(:), ALLOCATABLE :: arr_i4_1d
        INTEGER(8), DIMENSION(:), ALLOCATABLE :: arr_i8_1d
        INTEGER(1), DIMENSION(:), ALLOCATABLE :: arr_c1_1d
        REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: arr_r4_3d
        REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: arr_r8_3d
        REAL(16), DIMENSION(:,:,:), ALLOCATABLE :: arr_r16_3d
        INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE :: arr_i4_3d
        INTEGER(8), DIMENSION(:,:,:), ALLOCATABLE :: arr_i8_3d
        INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: arr_c1_3d

        inquire(file=B_filename, exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open B field file (from Lare3d)'
          print *, 'Attemped to read filename: '
          print *, TRIM(B_filename)
          call EXIT(113)
        END IF
        open(newunit=L3D_unit,file=B_filename,access='stream')
      !Read header
        read(L3D_unit) sdf1_str
        read(L3D_unit) endianness
        read(L3D_unit) sdf_version
        read(L3D_unit) sdf_revision
        read(L3D_unit) code_name
        read(L3D_unit) first_block_location
        read(L3D_unit) summary_location
        read(L3D_unit) summary_size
        read(L3D_unit) nblocks
        read(L3D_unit) block_header_length
        read(L3D_unit) step
        read(L3D_unit) time
        read(L3D_unit) jobid1
        read(L3D_unit) jobid2
        read(L3D_unit) string_length
        read(L3D_unit) code_io_version
        read(L3D_unit) restart_flag
        read(L3D_unit) subdomain_file
        ALLOCATE(character(string_length) :: block_name)

      !Read grid
        next_block_location = first_block_location
        do idx_b=1,nblocks
          call fseek(L3D_unit, next_block_location, 0)
          current_block_location = next_block_location
          read(L3D_unit) next_block_location
          read(L3D_unit) data_location
          read(L3D_unit) block_id
          if ((block_id(1:4) .eq. 'Grid') .or. (block_id(1:4) .eq. 'grid')) then
            read(L3D_unit) data_length
            read(L3D_unit) blocktype
            read(L3D_unit) datatype
            read(L3D_unit) ndims
            read(L3D_unit) block_name
            read(L3D_unit) block_info_length
            call fseek(L3D_unit, current_block_location+block_header_length, 0)
            ALLOCATE(mults(ndims))
            ALLOCATE(dims(ndims))
            read(L3D_unit) mults
            call fseek(L3D_unit, current_block_location+block_header_length+88*ndims+4, 0)
            read(L3D_unit) dims
            sz_1=dims(1)-1
            sz_2=dims(2)-1
            sz_3=dims(3)-1
            ALLOCATE(grid1(sz_1))
            ALLOCATE(grid2(sz_2))
            ALLOCATE(grid3(sz_3))
            ALLOCATE(B_grid(3,sz_1,sz_2,sz_3))
            call fseek(L3D_unit, data_location, 0)
            if (datatype .eq. 1) then
              ALLOCATE(arr_i4_1d(dims(1)))
              read(L3D_unit) arr_i4_1d
              grid1(:) = 0.5_num*(arr_i4_1d(1:dims(1)-1)+arr_i4_1d(2:dims(1)))
              DEALLOCATE(arr_i4_1d)
              ALLOCATE(arr_i4_1d(dims(2)))
              read(L3D_unit) arr_i4_1d
              grid2(:) = 0.5_num*(arr_i4_1d(1:dims(2)-1)+arr_i4_1d(2:dims(2)))
              DEALLOCATE(arr_i4_1d)
              ALLOCATE(arr_i4_1d(dims(3)))
              read(L3D_unit) arr_i4_1d
              grid3(:) = 0.5_num*(arr_i4_1d(1:dims(3)-1)+arr_i4_1d(2:dims(3)))
              DEALLOCATE(arr_i4_1d)
            else if (datatype .eq. 2) then
              ALLOCATE(arr_i8_1d(dims(1)))
              read(L3D_unit) arr_i8_1d
              grid1(:) = 0.5_num*(arr_i8_1d(1:dims(1)-1)+arr_i8_1d(2:dims(1)))
              DEALLOCATE(arr_i8_1d)
              ALLOCATE(arr_i8_1d(dims(2)))
              read(L3D_unit) arr_i8_1d
              grid2(:) = 0.5_num*(arr_i8_1d(1:dims(2)-1)+arr_i8_1d(2:dims(2)))
              DEALLOCATE(arr_i8_1d)
              ALLOCATE(arr_i8_1d(dims(3)))
              read(L3D_unit) arr_i8_1d
              grid3(:) = 0.5_num*(arr_i8_1d(1:dims(3)-1)+arr_i8_1d(2:dims(3)))
              DEALLOCATE(arr_i8_1d)
            else if (datatype .eq. 3) then
              ALLOCATE(arr_r4_1d(dims(1)))
              read(L3D_unit) arr_r4_1d
              grid1(:) = 0.5_num*(arr_r4_1d(1:dims(1)-1)+arr_r4_1d(2:dims(1)))
              DEALLOCATE(arr_r4_1d)
              ALLOCATE(arr_r4_1d(dims(2)))
              read(L3D_unit) arr_r4_1d
              grid2(:) = 0.5_num*(arr_r4_1d(1:dims(2)-1)+arr_r4_1d(2:dims(2)))
              DEALLOCATE(arr_r4_1d)
              ALLOCATE(arr_r4_1d(dims(3)))
              read(L3D_unit) arr_r4_1d
              grid3(:) = 0.5_num*(arr_r4_1d(1:dims(3)-1)+arr_r4_1d(2:dims(3)))
              DEALLOCATE(arr_r4_1d)
            else if (datatype .eq. 4) then
              ALLOCATE(arr_r8_1d(dims(1)))
              read(L3D_unit) arr_r8_1d
              grid1(:) = 0.5_num*(arr_r8_1d(1:dims(1)-1)+arr_r8_1d(2:dims(1)))
              DEALLOCATE(arr_r8_1d)
              ALLOCATE(arr_r8_1d(dims(2)))
              read(L3D_unit) arr_r8_1d
              grid2(:) = 0.5_num*(arr_r8_1d(1:dims(2)-1)+arr_r8_1d(2:dims(2)))
              DEALLOCATE(arr_r8_1d)
              ALLOCATE(arr_r8_1d(dims(3)))
              read(L3D_unit) arr_r8_1d
              grid3(:) = 0.5_num*(arr_r8_1d(1:dims(3)-1)+arr_r8_1d(2:dims(3)))
              DEALLOCATE(arr_r8_1d)
            else if (datatype .eq. 5) then
              ALLOCATE(arr_r16_1d(dims(1)))
              read(L3D_unit) arr_r16_1d
              grid1(:) = 0.5_num*(arr_r16_1d(1:dims(1)-1)+arr_r16_1d(2:dims(1)))
              DEALLOCATE(arr_r16_1d)
              ALLOCATE(arr_r16_1d(dims(2)))
              read(L3D_unit) arr_r16_1d
              grid2(:) = 0.5_num*(arr_r16_1d(1:dims(2)-1)+arr_r16_1d(2:dims(2)))
              DEALLOCATE(arr_r16_1d)
              ALLOCATE(arr_r16_1d(dims(3)))
              read(L3D_unit) arr_r16_1d
              grid3(:) = 0.5_num*(arr_r16_1d(1:dims(3)-1)+arr_r16_1d(2:dims(3)))
              DEALLOCATE(arr_r16_1d)
            else if ((datatype .eq. 6) .or. (datatype .eq. 7)) then
              ALLOCATE(arr_c1_1d(dims(1)))
              read(L3D_unit) arr_c1_1d
              grid1(:) = 0.5_num*(arr_c1_1d(1:dims(1)-1)+arr_c1_1d(2:dims(1)))
              DEALLOCATE(arr_c1_1d)
              ALLOCATE(arr_c1_1d(dims(2)))
              read(L3D_unit) arr_c1_1d
              grid2(:) = 0.5_num*(arr_c1_1d(1:dims(2)-1)+arr_c1_1d(2:dims(2)))
              DEALLOCATE(arr_c1_1d)
              ALLOCATE(arr_c1_1d(dims(3)))
              read(L3D_unit) arr_c1_1d
              grid3(:) = 0.5_num*(arr_c1_1d(1:dims(3)-1)+arr_c1_1d(2:dims(3)))
              DEALLOCATE(arr_c1_1d)
            end if
            DEALLOCATE(dims)
            DEALLOCATE(mults)
          end if
        end do

      !Read B field, regularize onto cell centers
        next_block_location = first_block_location
        do idx_b=1,nblocks
          call fseek(L3D_unit, next_block_location, 0)
          current_block_location = next_block_location
          read(L3D_unit) next_block_location
          read(L3D_unit) data_location
          read(L3D_unit) block_id
          if ((block_id(1:2) .eq. 'Bx') .or. (block_id(1:2) .eq. 'By') .or. &
                   (block_id(1:2) .eq. 'Bz')) then
            read(L3D_unit) data_length
            read(L3D_unit) blocktype
            read(L3D_unit) datatype
            read(L3D_unit) ndims
            read(L3D_unit) block_name
            read(L3D_unit) block_info_length
            call fseek(L3D_unit, current_block_location+block_header_length, 0)
            ALLOCATE(mults(ndims))
            ALLOCATE(dims(ndims))
            read(L3D_unit) mults
            call fseek(L3D_unit, current_block_location+block_header_length+72, 0)
            read(L3D_unit) dims
            call fseek(L3D_unit, data_location, 0)
            if (datatype .eq. 1) then
              ALLOCATE(arr_i4_3d(dims(1),dims(2),dims(3)))
              read(L3D_unit) arr_i4_3d
              if (block_id(1:2) .eq. 'Bx') then
                B_grid(1,:,:,:) = 0.5_num*(arr_i4_3d(1:dims(1)-1,:,:)+arr_i4_3d(2:dims(1),:,:))
              else if (block_id(1:2) .eq. 'By') then
                B_grid(2,:,:,:) = 0.5_num*(arr_i4_3d(:,1:dims(2)-1,:)+arr_i4_3d(:,2:dims(2),:))
              else if (block_id(1:2) .eq. 'Bz') then
                B_grid(3,:,:,:) = 0.5_num*(arr_i4_3d(:,:,1:dims(3)-1)+arr_i4_3d(:,:,2:dims(3)))
              end if
              DEALLOCATE(arr_i4_3d)
            else if (datatype .eq. 2) then
              ALLOCATE(arr_i8_3d(dims(1),dims(2),dims(3)))
              read(L3D_unit) arr_i8_3d
              if (block_id(1:2) .eq. 'Bx') then
                B_grid(1,:,:,:) = 0.5_num*(arr_i8_3d(1:dims(1)-1,:,:)+arr_i8_3d(2:dims(1),:,:))
              else if (block_id(1:2) .eq. 'By') then
                B_grid(2,:,:,:) = 0.5_num*(arr_i8_3d(:,1:dims(2)-1,:)+arr_i8_3d(:,2:dims(2),:))
              else if (block_id(1:2) .eq. 'Bz') then
                B_grid(3,:,:,:) = 0.5_num*(arr_i8_3d(:,:,1:dims(3)-1)+arr_i8_3d(:,:,2:dims(3)))
              end if
              DEALLOCATE(arr_i8_3d)
            else if (datatype .eq. 3) then
              ALLOCATE(arr_r4_3d(dims(1),dims(2),dims(3)))
              read(L3D_unit) arr_r4_3d
              if (block_id(1:2) .eq. 'Bx') then
                B_grid(1,:,:,:) = 0.5_num*(arr_r4_3d(1:dims(1)-1,:,:)+arr_r4_3d(2:dims(1),:,:))
              else if (block_id(1:2) .eq. 'By') then
                B_grid(2,:,:,:) = 0.5_num*(arr_r4_3d(:,1:dims(2)-1,:)+arr_r4_3d(:,2:dims(2),:))
              else if (block_id(1:2) .eq. 'Bz') then
                B_grid(3,:,:,:) = 0.5_num*(arr_r4_3d(:,:,1:dims(3)-1)+arr_r4_3d(:,:,2:dims(3)))
              end if
              DEALLOCATE(arr_r4_3d)
            else if (datatype .eq. 4) then
              ALLOCATE(arr_r8_3d(dims(1),dims(2),dims(3)))
              read(L3D_unit) arr_r8_3d
              if (block_id(1:2) .eq. 'Bx') then
                B_grid(1,:,:,:) = 0.5_num*(arr_r8_3d(1:dims(1)-1,:,:)+arr_r8_3d(2:dims(1),:,:))
              else if (block_id(1:2) .eq. 'By') then
                B_grid(2,:,:,:) = 0.5_num*(arr_r8_3d(:,1:dims(2)-1,:)+arr_r8_3d(:,2:dims(2),:))
              else if (block_id(1:2) .eq. 'Bz') then
                B_grid(3,:,:,:) = 0.5_num*(arr_r8_3d(:,:,1:dims(3)-1)+arr_r8_3d(:,:,2:dims(3)))
              end if
              DEALLOCATE(arr_r8_3d)
            else if (datatype .eq. 5) then
              ALLOCATE(arr_r16_3d(dims(1),dims(2),dims(3)))
              read(L3D_unit) arr_r16_3d
              if (block_id(1:2) .eq. 'Bx') then
                B_grid(1,:,:,:) = 0.5_num*(arr_r16_3d(1:dims(1)-1,:,:)+arr_r16_3d(2:dims(1),:,:))
              else if (block_id(1:2) .eq. 'By') then
                B_grid(2,:,:,:) = 0.5_num*(arr_r16_3d(:,1:dims(2)-1,:)+arr_r16_3d(:,2:dims(2),:))
              else if (block_id(1:2) .eq. 'Bz') then
                B_grid(3,:,:,:) = 0.5_num*(arr_r16_3d(:,:,1:dims(3)-1)+arr_r16_3d(:,:,2:dims(3)))
              end if
              DEALLOCATE(arr_r16_3d)
            else if ((datatype .eq. 6) .or. (datatype .eq. 7)) then
              ALLOCATE(arr_c1_3d(dims(1),dims(2),dims(3)))
              read(L3D_unit) arr_c1_3d
              if (block_id(1:2) .eq. 'Bx') then
                B_grid(1,:,:,:) = 0.5_num*(arr_c1_3d(1:dims(1)-1,:,:)+arr_c1_3d(2:dims(1),:,:))
              else if (block_id(1:2) .eq. 'By') then
                B_grid(2,:,:,:) = 0.5_num*(arr_c1_3d(:,1:dims(2)-1,:)+arr_c1_3d(:,2:dims(2),:))
              else if (block_id(1:2) .eq. 'Bz') then
                B_grid(3,:,:,:) = 0.5_num*(arr_c1_3d(:,:,1:dims(3)-1)+arr_c1_3d(:,:,2:dims(3)))
              end if
              DEALLOCATE(arr_c1_3d)
            end if
            DEALLOCATE(dims)
            DEALLOCATE(mults)
          end if
        end do

        close(L3D_unit)
        DEALLOCATE(block_name)

      end subroutine load_Lare3d


      subroutine load_ARMS_flicks

        INTEGER :: stp_idx, hdr_unit, flicks_unit, stat, stat2, qtynum, B_start_idx
        INTEGER :: idx_leaf, idx_blk, idx1, idx2, idx3, real_size
        LOGICAL :: log_grid, stop_found, B_field_found, bfile_exists
        INTEGER(4) :: ntblks, nlblks, newgrd, nvar
        INTEGER(4) :: iputwrk(21)
        REAL(4) :: time32
        REAL(4) :: rputwrk32(6)
        REAL(4), DIMENSION(:), ALLOCATABLE :: data_temp32
        REAL(8) :: time64
        REAL(8) :: rputwrk64(6)
        REAL(8), DIMENSION(:), ALLOCATABLE :: data_temp64
        CHARACTER(len=str_mx) :: hdr_line
        CHARACTER(len=4) :: blank4
        CHARACTER(len=8) :: blank8
        CHARACTER(len=25) :: blank25
        CHARACTER(len=str_mx) :: hdr_filename

        B_start_idx = 0
        stop_found = .false.
        B_field_found = .false.
        nvar = 0

        stp_idx = str_mx
        DO WHILE ((stp_idx .gt. 0) .and. (.not. stop_found))
          if (B_filename(stp_idx:stp_idx) .eq. '.') then
            stop_found = .true.
          else
            stp_idx = stp_idx - 1
          end if
        END DO
        hdr_filename=B_filename(1:stp_idx) // 'hdr'

        inquire(file=TRIM(hdr_filename), exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open flicks header file (must end with .hdr)'
          print *, 'This should be in the same directory as the target flicks file'
          print *, 'Attemped to read filename: '
          print *, TRIM(hdr_filename)
          call EXIT(130)
        END IF
        open(newunit=hdr_unit,file=TRIM(hdr_filename),form="formatted")
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        if (hdr_line(6:7) .eq. '32') then
          real_size = 4
        else if (hdr_line(6:7) .eq. '64') then
          real_size = 8
        else
          print *, 'Following line could not be interpreted as floating point size:'
          print *, hdr_line
          call EXIT(134)
        end if
        read(hdr_unit, '(A)', iostat=stat) hdr_line
      !Currently defined coordinate types
        if (TRIM(hdr_line) .eq. 'Cartesian') then
          geometry = 0
          log_grid = .false.
        else if (TRIM(hdr_line) .eq. 'Spherical Linear') then
          geometry = 1
          log_grid = .false.
        else if (TRIM(hdr_line) .eq. 'Spherical Exponential') then
          geometry = 1
          log_grid = .true.
        else if (TRIM(hdr_line) .eq. 'Cylindrical') then
          geometry = 2
          log_grid = .false.
        else
          print *, 'unable to identify ARMS coordinate type: ', TRIM(hdr_line)
          call EXIT(131)
        end if
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_line(1:2),*,iostat=stat2) sz_1
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_line(1:2),*,iostat=stat2) sz_2
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_line(1:2),*,iostat=stat2) sz_3
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        do while (stat .eq. 0)
          read(hdr_line(1:1),*,iostat=stat2) qtynum
          if (TRIM(hdr_line(2:str_mx)) .eq. 'Magnetic Field') then
            B_start_idx = nvar + 1
            B_field_found=.true.
          end if
          nvar = nvar + qtynum
          read(hdr_unit, '(A)', iostat=stat) hdr_line
        end do
        close(hdr_unit)
        if (.not. B_field_found) then
          print *, 'ARMS output does not contain Magnetic Field'
          call EXIT(132)
        end if
        if (real_size .eq. 4) then
          ALLOCATE(data_temp32(nvar))
        else if (real_size .eq. 8) then
          ALLOCATE(data_temp64(nvar))
        end if

        inquire(file=B_filename, exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open flicks file'
          print *, 'Attemped to read filename: '
          print *, TRIM(hdr_filename)
          call EXIT(133)
        END IF
        open(newunit=flicks_unit,file=B_filename,access='stream',convert='big_endian')
        read(flicks_unit) blank25
        if (real_size .eq. 4) then
          read(flicks_unit) time32
        else if (real_size .eq. 8) then
          read(flicks_unit) time64
        end if
        read(flicks_unit) blank8
        read(flicks_unit) ntblks
        read(flicks_unit) nlblks
        num_blocks = nlblks
        read(flicks_unit) newgrd
        read(flicks_unit) blank4
        ALLOCATE(grid1_ir(2,num_blocks))
        ALLOCATE(grid2_ir(2,num_blocks))
        ALLOCATE(grid3_ir(2,num_blocks))
        ALLOCATE(B_grid_ir(3,sz_1,sz_2,sz_3,num_blocks))

        if (real_size .eq. 4) then
          idx_leaf = 1
          DO idx_blk = 1,ntblks
            read(flicks_unit) blank4
            read(flicks_unit) iputwrk
            read(flicks_unit) blank8
            read(flicks_unit) rputwrk32
            read(flicks_unit) blank4
            if (iputwrk(3) .eq. 1) then
              if (log_grid) then
                grid1_ir(:,idx_leaf) = EXP(rputwrk32(1:2))
              else
                grid1_ir(:,idx_leaf) = rputwrk32(1:2)
              end if
              if (geometry .eq. 1) then
                grid2_ir(:,idx_leaf) = 0.5_num*PI - rputwrk32(3:4)
              else
                grid2_ir(:,idx_leaf) = rputwrk32(3:4)
              end if
              if (geometry .eq. 1) then
                grid3_ir(:,idx_leaf) = rputwrk32(5:6) + PI
              else
                grid3_ir(:,idx_leaf) = rputwrk32(5:6)
              end if
              DO idx3 = 1,sz_3
                DO idx2 = 1,sz_2
                  DO idx1 = 1,sz_1
                    read(flicks_unit) blank4
                    read(flicks_unit) data_temp32
                    B_grid_ir(:,idx1,idx2,idx3,idx_leaf)=data_temp32(B_start_idx:B_start_idx+2)
                    read(flicks_unit) blank4
                  END DO
                END DO
              END DO
              idx_leaf = idx_leaf + 1
            end if
          END DO
        else if (real_size .eq. 8) then
          idx_leaf = 1
          DO idx_blk = 1,ntblks
            read(flicks_unit) blank4
            read(flicks_unit) iputwrk
            read(flicks_unit) blank8
            read(flicks_unit) rputwrk64
            read(flicks_unit) blank4
            if (iputwrk(3) .eq. 1) then
              if (log_grid) then
                grid1_ir(:,idx_leaf) = EXP(rputwrk64(1:2))
              else
                grid1_ir(:,idx_leaf) = rputwrk64(1:2)
              end if
              if (geometry .eq. 1) then
                grid2_ir(:,idx_leaf) = 0.5_num*PI - rputwrk64(3:4)
              else
                grid2_ir(:,idx_leaf) = rputwrk64(3:4)
              end if
              if (geometry .eq. 1) then
                grid3_ir(:,idx_leaf) = rputwrk64(5:6) + PI
              else
                grid3_ir(:,idx_leaf) = rputwrk64(5:6)
              end if
              DO idx3 = 1,sz_3
                DO idx2 = 1,sz_2
                  DO idx1 = 1,sz_1
                    read(flicks_unit) blank4
                    read(flicks_unit) data_temp64
                    B_grid_ir(:,idx1,idx2,idx3,idx_leaf)=data_temp64(B_start_idx:B_start_idx+2)
                    read(flicks_unit) blank4
                  END DO
                END DO
              END DO
              idx_leaf = idx_leaf + 1
            end if
          END DO
        end if

        close(flicks_unit)
        if (real_size .eq. 4) then
          DEALLOCATE(data_temp32)
        else if (real_size .eq. 8) then
          DEALLOCATE(data_temp64)
        end if

      end subroutine load_ARMS_flicks


      subroutine load_ARMS_bfield

        INTEGER :: stp_idx, hdr_unit, bfield_unit, stat, stat2
        INTEGER :: idx_leaf, idx_blk, idx1, idx2, idx3
        LOGICAL :: log_grid, stop_found, bfile_exists
        INTEGER(4) :: ntblks, nlblks
        INTEGER(4) :: iputwrk(35)
        REAL(8) :: time64
        REAL(8) :: rputwrk64(6)
        REAL(8), DIMENSION(:), ALLOCATABLE :: data_temp64
        CHARACTER(len=str_mx) :: hdr_line
        CHARACTER(len=4) :: blank4
        CHARACTER(len=8) :: blank8
        CHARACTER(len=str_mx) :: hdr_filename

        stop_found = .false.

        stp_idx = str_mx
        DO WHILE ((stp_idx .gt. 0) .and. (.not. stop_found))
          if (B_filename(stp_idx:stp_idx) .eq. '.') then
            stop_found = .true.
          else
            stp_idx = stp_idx - 1
          end if
        END DO
        hdr_filename=B_filename(1:stp_idx) // 'hdr'

        inquire(file=TRIM(hdr_filename), exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open bfield header file (must end with .hdr)'
          print *, 'This should be in the same directory as the target bfield file'
          print *, 'Attemped to read filename: '
          print *, TRIM(hdr_filename)
          call EXIT(135)
        END IF
        open(newunit=hdr_unit,file=TRIM(hdr_filename),form="formatted")
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_unit, '(A)', iostat=stat) hdr_line
      !Currently defined coordinate types
        if (TRIM(hdr_line) .eq. 'Cartesian') then
          geometry = 0
          log_grid = .false.
        else if (TRIM(hdr_line) .eq. 'Spherical Linear') then
          geometry = 1
          log_grid = .false.
        else if (TRIM(hdr_line) .eq. 'Spherical Exponential') then
          geometry = 1
          log_grid = .true.
        else if (TRIM(hdr_line) .eq. 'Cylindrical') then
          geometry = 2
          log_grid = .false.
        else
          print *, 'unable to identify ARMS coordinate type: ', TRIM(hdr_line)
          call EXIT(136)
        end if
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_line(1:2),*,iostat=stat2) sz_1
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_line(1:2),*,iostat=stat2) sz_2
        read(hdr_unit, '(A)', iostat=stat) hdr_line
        read(hdr_line(1:2),*,iostat=stat2) sz_3
        close(hdr_unit)

        ALLOCATE(data_temp64(3))
        inquire(file=B_filename, exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open bfield file'
          print *, 'Attemped to read filename: '
          print *, TRIM(hdr_filename)
          call EXIT(137)
        END IF
        open(newunit=bfield_unit,file=B_filename,access='stream',convert='big_endian')
        read(bfield_unit) blank4
        read(bfield_unit) time64
        read(bfield_unit) blank8
        read(bfield_unit) ntblks
        read(bfield_unit) nlblks
        num_blocks = nlblks
        read(bfield_unit) blank4
        ALLOCATE(grid1_ir(2,num_blocks))
        ALLOCATE(grid2_ir(2,num_blocks))
        ALLOCATE(grid3_ir(2,num_blocks))
        ALLOCATE(B_grid_ir(3,sz_1,sz_2,sz_3,num_blocks))

        idx_leaf = 1
        DO idx_blk = 1,ntblks
          read(bfield_unit) blank4
          read(bfield_unit) iputwrk
          read(bfield_unit) blank8
          read(bfield_unit) rputwrk64
          read(bfield_unit) blank4
          if (iputwrk(3) .eq. 1) then
            if (log_grid) then
              grid1_ir(:,idx_leaf) = EXP(rputwrk64(1:2))
            else
              grid1_ir(:,idx_leaf) = rputwrk64(1:2)
            end if
            if (geometry .eq. 1) then
              grid2_ir(:,idx_leaf) = 0.5_num*PI - rputwrk64(3:4)
            else
              grid2_ir(:,idx_leaf) = rputwrk64(3:4)
            end if
            if (geometry .eq. 1) then
              grid3_ir(:,idx_leaf) = rputwrk64(5:6) + PI
            else
              grid3_ir(:,idx_leaf) = rputwrk64(5:6)
            end if
            DO idx3 = 1,sz_3
              DO idx2 = 1,sz_2
                DO idx1 = 1,sz_1
                  read(bfield_unit) blank4
                  read(bfield_unit) data_temp64
                  B_grid_ir(:,idx1,idx2,idx3,idx_leaf)=data_temp64
                  read(bfield_unit) blank4
                END DO
              END DO
            END DO
            idx_leaf = idx_leaf + 1
          end if
        END DO

        close(bfield_unit)
        DEALLOCATE(data_temp64)

      end subroutine load_ARMS_bfield


      subroutine write_output

        INTEGER :: out_unit, idx
        CHARACTER(HEADER_LENGTH) :: header_start
        CHARACTER(HEADER_LENGTH) :: string_temp

        header_start = 'UFiT Output, v1.0, Header length: '
        write(string_temp, '(I0)') HEADER_LENGTH
        header_start = trim(header_start) // ' ' // trim(string_temp) // ';'
        write(string_temp, '(I0)') num
        header_start = trim(header_start) // ' N: ' // trim(string_temp) // ';'
        write(string_temp, '(I0)') geometry
        header_start = trim(header_start) // ' G: ' // trim(string_temp) // ';'
        write(string_temp, '(I0)') input_type
        header_start = trim(header_start) // ' IT: ' // trim(string_temp) // ';'
        write(string_temp, *) include_curvature
        header_start = trim(header_start) // ' IC:' // trim(string_temp) // ';'
        write(string_temp, *) save_endpoints
        header_start = trim(header_start) // ' SE:' // trim(string_temp) // ';'
        write(string_temp, *) save_Q
        header_start = trim(header_start) // ' SQ:' // trim(string_temp) // ';'
        write(string_temp, *) save_fieldlines
        header_start = trim(header_start) // ' SF:' // trim(string_temp) // ';'
        write(string_temp, *) save_connection
        header_start = trim(header_start) // ' SC:' // trim(string_temp) // ';'
        write(string_temp, *) user_defined
        header_start = trim(header_start) // ' UD:' // trim(string_temp) // ';'
        write(string_temp, *) periodic_X
        header_start = trim(header_start) // ' PX:' // trim(string_temp) // ';'
        write(string_temp, *) periodic_Y
        header_start = trim(header_start) // ' PY:' // trim(string_temp) // ';'
        write(string_temp, *) periodic_Z
        header_start = trim(header_start) // ' PZ:' // trim(string_temp) // ';'
        write(string_temp, *) periodic_PHI
        header_start = trim(header_start) // ' PP:' // trim(string_temp) // ';'
        write(string_temp, *) grid_separate
        header_start = trim(header_start) // ' GS:' // trim(string_temp) // ';'

        open(newunit=out_unit, file=out_filename ,access = 'stream', status = 'replace')
        write(out_unit) header_start
      !Write input grid
        if (input_type .eq. 0) then
          write(out_unit) numin_tot
          write(out_unit) coord1_in
          write(out_unit) coord2_in
          write(out_unit) coord3_in
        else !input_type == 1 or input_type == 2
          write(out_unit) numin1
          write(out_unit) numin2
          write(out_unit) numin3
          write(out_unit) coord1_in
          write(out_unit) coord2_in
          write(out_unit) coord3_in
        end if

      !Write endpoints
        if (save_endpoints) then
          do idx = 1, numin_tot
            write(out_unit) fieldline_endpoints(:,idx)
          end do
        end if
      !Write Q
        if (save_Q) then
          do idx = 1, numin_tot
            write(out_unit) fieldline_Q(idx)
          end do
        end if
      !Write fieldlines
        if (save_fieldlines) then
          do idx = 1, numin_tot
            write(out_unit) fieldline_ptn(idx)
          end do
          do idx = 1, numin_tot
            write(out_unit) fieldline_allpos(:,fieldline_pts(idx):fieldline_pts(idx) &
                                             +fieldline_ptn(idx)-1,idx)
          end do
        end if
      !Write connection
        if (save_connection) then
          do idx = 1, numin_tot
            write(out_unit) fieldline_connection(idx)
          end do
        end if
      !User defined
        if (user_defined) then
          write(out_unit) num_ud_variables
          do idx = 1, numin_tot
            write(out_unit) fieldline_user(:,idx)
          end do
        end if

        close(out_unit)

      end subroutine write_output


      subroutine get_start_pos_0(pos_start,idx_in)
      !Input type 0 = each point is explicitly by set of coordinates

        REAL(num) :: pos_start(3)
        INTEGER :: idx_in

        pos_start(1) = coord1_in(idx_in)
        pos_start(2) = coord2_in(idx_in)
        pos_start(3) = coord3_in(idx_in)

      end subroutine get_start_pos_0


      subroutine get_start_pos_12(pos_start,idx_in)
      !Input type 1 = seed coordinates are a regular grid; each axis is specified explicitly
      !        or 2 = seed coordinates are an evenly-spaced grid, specified by the start and end

        REAL(num) :: pos_start(3)
        INTEGER :: idx_in

        INTEGER :: idx1, idx2, idx3

        idx1=MOD(idx_in-1,numin1)+1
        idx2=MOD((idx_in-1)/numin1,numin2)+1
        idx3=(idx_in-1)/(numin2*numin1)+1

        pos_start(1) = coord1_in(idx1)
        pos_start(2) = coord2_in(idx2)
        pos_start(3) = coord3_in(idx3)

      end subroutine get_start_pos_12


      function find_index(grid,sz,target_loc,idx_in)
      !Index of gridpoint below ; grid is increasing

        INTEGER :: sz, idx_in
        REAL(num) :: target_loc
        REAL(num) :: grid(sz)

        INTEGER :: find_index

        find_index = idx_in
        if (grid(find_index) .le. target_loc) then
          do while ((find_index .lt. sz) .and. (grid(find_index+1) .le. target_loc))
            find_index = find_index + 1
          end do
        else
          do while ((find_index .gt. 1) .and. (grid(find_index) .gt. target_loc))
            find_index = find_index - 1
          end do
        end if

        find_index=MAX(MIN(find_index,sz-1),1)

      end function find_index


      function find_index_irregular(target_loc,idx_in)
      !Index of gridpoint below ; grid is increasing

        INTEGER :: idx_in
        REAL(num) :: target_loc(3)

        INTEGER :: idx
        INTEGER :: find_index_irregular

        if ((target_loc(1) .ge. grid1_ir(1,idx_in)) .and. (target_loc(1) .lt. grid1_ir(2,idx_in)) &
           .and. (target_loc(2) .ge. grid2_ir(1,idx_in)) .and. (target_loc(2) .lt. &
           grid2_ir(2,idx_in)) .and. (target_loc(3) .ge. grid3_ir(1,idx_in)) .and. (target_loc(3) &
           .lt. grid3_ir(2,idx_in))) then
          find_index_irregular = idx_in
        else
          do idx=1,num_blocks
            if ((target_loc(1) .ge. grid1_ir(1,idx)) .and. (target_loc(1) .lt. &
                grid1_ir(2,idx)) .and. (target_loc(2) .ge. grid2_ir(1,idx)) .and. &
                (target_loc(2) .lt. grid2_ir(2,idx)) .and. (target_loc(3) .ge. &
                grid3_ir(1,idx)) .and. (target_loc(3) .lt. grid3_ir(2,idx))) then
              find_index_irregular = idx
              exit
            end if
          end do
        end if

      end function find_index_irregular


      function vecdot(v1,v2)

        REAL(num) :: v1(3), v2(3)

        REAL(num) :: vecdot

        vecdot = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)

      end function vecdot


      function normalize_vector(vec_in)

        REAL(num) :: vec_in(3)

        REAL(num) :: normalize_vector(3)

        normalize_vector(:) = vec_in(:)/SQRT(vec_in(1)**2+vec_in(2)**2+vec_in(3)**2)

      end function normalize_vector


      subroutine check_position_c000(pos_in,pos_out,keep_running)
      !Cartesian, non-periodic in X, Y, Z

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        if (pos_in(3) .lt. grid3min) then
          pos_out(3) = grid3min
          keep_running = .false.
        else if (pos_in(3) .gt. grid3max) then
          pos_out(3) = grid3max
          keep_running = .false.
        end if

      end subroutine check_position_c000


      subroutine check_position_c100(pos_in,pos_out,keep_running)
      !Cartesian, periodic in X, non-periodic in Y, Z

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(:) = pos_in(:)

        pos_out(1) = MODULO(pos_in(1)-grid1min,coord_width(1))+grid1min

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        if (pos_in(3) .lt. grid3min) then
          pos_out(3) = grid3min
          keep_running = .false.
        else if (pos_in(3) .gt. grid3max) then
          pos_out(3) = grid3max
          keep_running = .false.
        end if

      end subroutine check_position_c100


      subroutine check_position_c010(pos_in,pos_out,keep_running)
      !Cartesian, periodic in Y, non-periodic in X, Z

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        pos_out(2) = MODULO(pos_in(2)-grid2min,coord_width(2))+grid2min

        if (pos_in(3) .lt. grid3min) then
          pos_out(3) = grid3min
          keep_running = .false.
        else if (pos_in(3) .gt. grid3max) then
          pos_out(3) = grid3max
          keep_running = .false.
        end if

      end subroutine check_position_c010


      subroutine check_position_c001(pos_in,pos_out,keep_running)
      !Cartesian, periodic in Z, non-periodic in X, Y

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        pos_out(3) = MODULO(pos_in(3)-grid3min,coord_width(3))+grid3min

      end subroutine check_position_c001


      subroutine check_position_c110(pos_in,pos_out,keep_running)
      !Cartesian, periodic in X, Y non-periodic in Z

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(3) = pos_in(3)

        pos_out(1) = MODULO(pos_in(1)-grid1min,coord_width(1))+grid1min

        pos_out(2) = MODULO(pos_in(2)-grid2min,coord_width(2))+grid2min

        if (pos_in(3) .lt. grid3min) then
          pos_out(3) = grid3min
          keep_running = .false.
        else if (pos_in(3) .gt. grid3max) then
          pos_out(3) = grid3max
          keep_running = .false.
        end if

      end subroutine check_position_c110


      subroutine check_position_c101(pos_in,pos_out,keep_running)
      !Cartesian, periodic in X, Z, non-periodic in Y

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(2) = pos_in(2)

        pos_out(1) = MODULO(pos_in(1)-grid1min,coord_width(1))+grid1min

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        pos_out(3) = MODULO(pos_in(3)-grid3min,coord_width(3))+grid3min

      end subroutine check_position_c101


      subroutine check_position_c011(pos_in,pos_out,keep_running)
      !Cartesian, periodic in Y, Z, non-periodic in X

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(1) = pos_in(1)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        pos_out(2) = MODULO(pos_in(2)-grid2min,coord_width(2))+grid2min

        pos_out(3) = MODULO(pos_in(3)-grid3min,coord_width(3))+grid3min

      end subroutine check_position_c011


      subroutine check_position_c111(pos_in,pos_out,keep_running)
      !Cartesian, periodic in X, Y, Z

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(1) = MODULO(pos_in(1)-grid1min,coord_width(1))+grid1min

        pos_out(2) = MODULO(pos_in(2)-grid2min,coord_width(2))+grid2min

        pos_out(3) = MODULO(pos_in(3)-grid3min,coord_width(3))+grid3min

      end subroutine check_position_c111


      subroutine intercept_boundary_c000(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Cartesian, non-periodic in X, Y, Z

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(4)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(2)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(2) .lt. 0.0_num) then
          delta(2) = (grid2min - pos_in(2))/direction(2)
        else 
          delta(2) = (grid2max - pos_in(2))/direction(2)
        end if

        if (ABS(direction(3)) .lt. EPS) then
          delta(3) = 1.0_num/EPS
        else if (direction(3) .lt. 0.0_num) then
          delta(3) = (grid3min - pos_in(3))/direction(3)
        else 
          delta(3) = (grid3max - pos_in(3))/direction(3)
        end if

        delta(4) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_c000


      subroutine intercept_boundary_c100(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Cartesian, periodic in X, non-periodic in Y, Z

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(3)

        if (ABS(direction(2)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(2) .lt. 0.0_num) then
          delta(1) = (grid2min - pos_in(2))/direction(2)
        else 
          delta(1) = (grid2max - pos_in(2))/direction(2)
        end if

        if (ABS(direction(3)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(3) .lt. 0.0_num) then
          delta(2) = (grid3min - pos_in(3))/direction(3)
        else 
          delta(2) = (grid3max - pos_in(3))/direction(3)
        end if

        delta(3) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_c100


      subroutine intercept_boundary_c010(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Cartesian, periodic in Y, non-periodic in X, Z

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(3)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(3)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(3) .lt. 0.0_num) then
          delta(2) = (grid3min - pos_in(3))/direction(3)
        else 
          delta(2) = (grid3max - pos_in(3))/direction(3)
        end if

        delta(3) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_c010


      subroutine intercept_boundary_c001(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Cartesian, periodic in Z, non-periodic in X, Y

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(3)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(2)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(2) .lt. 0.0_num) then
          delta(2) = (grid2min - pos_in(2))/direction(2)
        else 
          delta(2) = (grid2max - pos_in(2))/direction(2)
        end if

        delta(3) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_c001


      subroutine intercept_boundary_c110(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Cartesian, periodic in X and Y, non-periodic in Z

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(2)

        if (ABS(direction(3)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(3) .lt. 0.0_num) then
          delta(1) = (grid3min - pos_in(3))/direction(3)
        else 
          delta(1) = (grid3max - pos_in(3))/direction(3)
        end if

        delta(2) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_c110


      subroutine intercept_boundary_c101(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Cartesian, periodic in X, Z, non-periodic in Y

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(2)

        if (ABS(direction(2)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(2) .lt. 0.0_num) then
          delta(1) = (grid2min - pos_in(2))/direction(2)
        else 
          delta(1) = (grid2max - pos_in(2))/direction(2)
        end if

        delta(2) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_c101


      subroutine intercept_boundary_c011(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Cartesian, periodic in Y, Z, non-periodic in X

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(2)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        delta(2) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_c011


      subroutine intercept_boundary_c111(pos_in,direction,delta_out)
      !Get distance to nearest boundary (doesn't exist in this case)
      !Cartesian, periodic in X, Y, Z

        REAL(num) :: pos_in(3),direction(3),delta_out

        delta_out = 1.0_num

      end subroutine intercept_boundary_c111


      subroutine check_position_s000(pos_in,pos_out,keep_running)
      !Spherical, non-periodic in r, theta, phi

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        if (pos_in(3) .lt. grid3min) then
          pos_out(3) = grid3min
          keep_running = .false.
        else if (pos_in(3) .gt. grid3max) then
          pos_out(3) = grid3max
          keep_running = .false.
        end if

      end subroutine check_position_s000


      subroutine check_position_s010(pos_in,pos_out,keep_running)
      !Spherical, non-periodic in r, phi
      !North theta boundary allows passing through

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        if (pos_in(2) .lt. 0.0_num) then
          pos_in(2) = -pos_in(2)
          if (pos_in(3) .lt. PI) then
            pos_in(3) = pos_in(3)+PI
          else
            pos_in(3) = pos_in(3)-PI
          end if
        else if (pos_in(2) .lt. grid2min) then
          pos_in(2) = grid2min
        end if
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        if (pos_in(3) .lt. grid3min) then
          pos_out(3) = grid3min
          keep_running = .false.
        else if (pos_in(3) .gt. grid3max) then
          pos_out(3) = grid3max
          keep_running = .false.
        end if

      end subroutine check_position_s010


      subroutine check_position_s020(pos_in,pos_out,keep_running)
      !Spherical, non-periodic in r, phi
      !South theta boundary allows passing through

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        if (pos_in(2) .gt. PI) then
          pos_in(2) = TWOPI-pos_in(2)
          if (pos_in(3) .lt. PI) then
            pos_in(3) = pos_in(3)+PI
          else
            pos_in(3) = pos_in(3)-PI
          end if
        else if (pos_in(2) .gt. grid2max) then
          pos_in(2) = grid2max
        end if
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        if (pos_in(3) .lt. grid3min) then
          pos_out(3) = grid3min
          keep_running = .false.
        else if (pos_in(3) .gt. grid3max) then
          pos_out(3) = grid3max
          keep_running = .false.
        end if

      end subroutine check_position_s020


      subroutine check_position_s030(pos_in,pos_out,keep_running)
      !Spherical, non-periodic in r, phi
      !Both theta boundaries allows passing through

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        if (pos_in(2) .lt. 0.0_num) then
          pos_in(2) = -pos_in(2)
          if (pos_in(3) .lt. PI) then
            pos_in(3) = pos_in(3)+PI
          else
            pos_in(3) = pos_in(3)-PI
          end if
        else if (pos_in(2) .lt. grid2min) then
          pos_in(2) = grid2min
        else if (pos_in(2) .gt. PI) then
          pos_in(2) = TWOPI-pos_in(2)
          if (pos_in(3) .lt. PI) then
            pos_in(3) = pos_in(3)+PI
          else
            pos_in(3) = pos_in(3)-PI
          end if
        else if (pos_in(2) .gt. grid2max) then
          pos_in(2) = grid2max
        end if
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        if (pos_in(3) .lt. grid3min) then
          pos_out(3) = grid3min
          keep_running = .false.
        else if (pos_in(3) .gt. grid3max) then
          pos_out(3) = grid3max
          keep_running = .false.
        end if

      end subroutine check_position_s030


      subroutine check_position_s001(pos_in,pos_out,keep_running)
      !Spherical, periodic in phi, non-periodic in r, theta

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        pos_out(3) = MODULO(pos_in(3)-grid3min,coord_width(3))+grid3min

      end subroutine check_position_s001


      subroutine check_position_s011(pos_in,pos_out,keep_running)
      !Spherical, periodic in phi, non-periodic in r
      !North theta boundary allows passing through

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        if (pos_in(2) .lt. 0.0_num) then
          pos_in(2) = -pos_in(2)
          if (pos_in(3) .lt. PI) then
            pos_in(3) = pos_in(3)+PI
          else
            pos_in(3) = pos_in(3)-PI
          end if
        else if (pos_in(2) .lt. grid2min) then
          pos_in(2) = grid2min
        end if
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        pos_out(3) = MODULO(pos_in(3)-grid3min,coord_width(3))+grid3min

      end subroutine check_position_s011


      subroutine check_position_s021(pos_in,pos_out,keep_running)
      !Spherical, periodic in phi, non-periodic in r
      !South theta boundary allows passing through

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        if (pos_in(2) .gt. PI) then
          pos_in(2) = TWOPI-pos_in(2)
          if (pos_in(3) .lt. PI) then
            pos_in(3) = pos_in(3)+PI
          else
            pos_in(3) = pos_in(3)-PI
          end if
        else if (pos_in(2) .gt. grid2max) then
          pos_in(2) = grid2max
        end if
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        pos_out(3) = MODULO(pos_in(3)-grid3min,coord_width(3))+grid3min

      end subroutine check_position_s021


      subroutine check_position_s031(pos_in,pos_out,keep_running)
      !Spherical, periodic in phi, non-periodic in r
      !Both theta boundaries allows passing through

        REAL(num) :: pos_in(3),pos_out(3)
        LOGICAL :: keep_running

        keep_running = .true.
        if (pos_in(2) .lt. 0.0_num) then
          pos_in(2) = -pos_in(2)
          if (pos_in(3) .lt. PI) then
            pos_in(3) = pos_in(3)+PI
          else
            pos_in(3) = pos_in(3)-PI
          end if
        else if (pos_in(2) .lt. grid2min) then
          pos_in(2) = grid2min
        else if (pos_in(2) .gt. PI) then
          pos_in(2) = TWOPI-pos_in(2)
          if (pos_in(3) .lt. PI) then
            pos_in(3) = pos_in(3)+PI
          else
            pos_in(3) = pos_in(3)-PI
          end if
        else if (pos_in(2) .gt. grid2max) then
          pos_in(2) = grid2max
        end if
        pos_out(:) = pos_in(:)

        if (pos_in(1) .lt. grid1min) then
          pos_out(1) = grid1min
          keep_running = .false.
        else if (pos_in(1) .gt. grid1max) then
          pos_out(1) = grid1max
          keep_running = .false.
        end if

        if (pos_in(2) .lt. grid2min) then
          pos_out(2) = grid2min
          keep_running = .false.
        else if (pos_in(2) .gt. grid2max) then
          pos_out(2) = grid2max
          keep_running = .false.
        end if

        pos_out(3) = MODULO(pos_in(3)-grid3min,coord_width(3))+grid3min

      end subroutine check_position_s031


      subroutine intercept_boundary_s000(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Spherical, non-periodic in r, theta, phi

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(4)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(2)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(2) .lt. 0.0_num) then
          delta(2) = (grid2min - pos_in(2))*pos_in(1)/direction(2)
        else 
          delta(2) = (grid2max - pos_in(2))*pos_in(1)/direction(2)
        end if

        if (ABS(direction(3)) .lt. EPS) then
          delta(3) = 1.0_num/EPS
        else if (direction(3) .lt. 0.0_num) then
          delta(3) = (grid3min - pos_in(3))*pos_in(1)*SIN(pos_in(2))/direction(3)
        else 
          delta(3) = (grid3max - pos_in(3))*pos_in(1)*SIN(pos_in(2))/direction(3)
        end if

        delta(4) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_s000


      subroutine intercept_boundary_s010(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Spherical, non-periodic in r, phi, theta open in North

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(4)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(2)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(2) .gt. 0.0_num) then
          delta(2) = (grid2max - pos_in(2))*pos_in(1)/direction(2)
        end if

        if (ABS(direction(3)) .lt. EPS) then
          delta(3) = 1.0_num/EPS
        else if (direction(3) .lt. 0.0_num) then
          delta(3) = (grid3min - pos_in(3))*pos_in(1)*SIN(pos_in(2))/direction(3)
        else 
          delta(3) = (grid3max - pos_in(3))*pos_in(1)*SIN(pos_in(2))/direction(3)
        end if

        delta(4) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_s010


      subroutine intercept_boundary_s020(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Spherical, non-periodic in r, phi, theta open in South

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(4)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(2)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(2) .lt. 0.0_num) then
          delta(2) = (grid2min - pos_in(2))*pos_in(1)/direction(2)
        end if

        if (ABS(direction(3)) .lt. EPS) then
          delta(3) = 1.0_num/EPS
        else if (direction(3) .lt. 0.0_num) then
          delta(3) = (grid3min - pos_in(3))*pos_in(1)*SIN(pos_in(2))/direction(3)
        else 
          delta(3) = (grid3max - pos_in(3))*pos_in(1)*SIN(pos_in(2))/direction(3)
        end if

        delta(4) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_s020


      subroutine intercept_boundary_s030(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Spherical, non-periodic in r, phi, theta open in both boundaries

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(3)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(3)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(3) .lt. 0.0_num) then
          delta(2) = (grid3min - pos_in(3))*pos_in(1)*SIN(pos_in(2))/direction(3)
        else 
          delta(2) = (grid3max - pos_in(3))*pos_in(1)*SIN(pos_in(2))/direction(3)
        end if

        delta(3) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_s030


      subroutine intercept_boundary_s001(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Spherical, periodic in phi, non-periodic in r, theta

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(3)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(2)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(2) .lt. 0.0_num) then
          delta(2) = (grid2min - pos_in(2))*pos_in(1)/direction(2)
        else 
          delta(2) = (grid2max - pos_in(2))*pos_in(1)/direction(2)
        end if

        delta(3) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_s001


      subroutine intercept_boundary_s011(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Spherical, non-periodic in r, phi, theta open in North

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(3)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(2)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(2) .gt. 0.0_num) then
          delta(2) = (grid2max - pos_in(2))*pos_in(1)/direction(2)
        end if

        delta(3) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_s011


      subroutine intercept_boundary_s021(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Spherical, non-periodic in r, phi, theta open in South

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(3)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        if (ABS(direction(2)) .lt. EPS) then
          delta(2) = 1.0_num/EPS
        else if (direction(2) .lt. 0.0_num) then
          delta(2) = (grid2min - pos_in(2))*pos_in(1)/direction(2)
        end if

        delta(3) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_s021


      subroutine intercept_boundary_s031(pos_in,direction,delta_out)
      !Get distance to nearest boundary
      !Spherical, non-periodic in r, phi, theta open in both boundaries

        REAL(num) :: pos_in(3),direction(3),delta_out

        REAL(num) :: delta(2)

        if (ABS(direction(1)) .lt. EPS) then
          delta(1) = 1.0_num/EPS
        else if (direction(1) .lt. 0.0_num) then
          delta(1) = (grid1min - pos_in(1))/direction(1)
        else 
          delta(1) = (grid1max - pos_in(1))/direction(1)
        end if

        delta(2) = 1.0_num

        delta_out = MINVAL(delta)

      end subroutine intercept_boundary_s031


      subroutine create_perpendicular_vectors(vec_in,u_out,v_out)
      !Create two normalized vectors perpendicular to input vector

        REAL(num) :: vec_in(3),u_out(3),v_out(3)
        REAL(num) :: vec_norm(3), vec_mod

        vec_norm(:) = vec_in(:)/SQRT(vec_in(1)**2+vec_in(2)**2+vec_in(3)**2)
        IF (vec_norm(3) .lt. EPS) THEN          !Effectively z-component is zero
          IF (vec_norm(1) .lt. EPS) THEN        !x-component also vanishingly small
            u_out(1) = 1.0_num
            u_out(2) = -vec_norm(1)/vec_norm(2)
            u_out(3) = 0.0_num
          ELSE                                  !x-component not too small
            u_out(1) = -vec_norm(2)/vec_norm(1)
            u_out(2) = 1.0_num
            u_out(3) = 0.0_num
          END IF
        ELSE   !z-component not too small
          u_out(1) = 0.0_num
          u_out(2) = 1.0_num
          u_out(3) = -vec_norm(2)/vec_norm(3)
        END IF

        !vec_in x u_out
        v_out(1) = vec_norm(2)*u_out(3)-vec_norm(3)*u_out(2)
        v_out(2) = vec_norm(3)*u_out(1)-vec_norm(1)*u_out(3)
        v_out(3) = vec_norm(1)*u_out(2)-vec_norm(2)*u_out(1)

        vec_mod = SQRT(u_out(1)**2+u_out(2)**2+u_out(3)**2)
        u_out = u_out(:) / vec_mod
        vec_mod = SQRT(v_out(1)**2+v_out(2)**2+v_out(3)**2)
        v_out = v_out(:) / vec_mod

      end subroutine create_perpendicular_vectors


      subroutine B_interp_regular(idx_in, pos_in, B_out)
      !Trilinear interpolation of B at a coordinate surrounded by grid cells
      !First transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)

        INTEGER :: idx1, idx2, idx3
        INTEGER :: idx1p, idx2p, idx3p
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1(3), a2(3), a3(3), a4(3), a5(3), a6(3), a7(3)
        REAL(num) :: B000(3), B001(3), B010(3), B011(3), B100(3), B101(3), B110(3), B111(3)

        idx1 = idx_in(1)
        idx2 = idx_in(2)
        idx3 = idx_in(3)
        idx1 = find_index(grid1,sz_1,pos_in(1),idx1)
        idx2 = find_index(grid2,sz_2,pos_in(2),idx2)
        idx3 = find_index(grid3,sz_3,pos_in(3),idx3)
        idx_in(1) = idx1
        idx_in(2) = idx2
        idx_in(3) = idx3
        idx1p = idx1+1
        idx2p = idx2+1
        idx3p = idx3+1
        B000(:) = B_grid(:,idx1,idx2,idx3)
        B001(:) = B_grid(:,idx1,idx2,idx3p)
        B010(:) = B_grid(:,idx1,idx2p,idx3)
        B011(:) = B_grid(:,idx1,idx2p,idx3p)
        B100(:) = B_grid(:,idx1p,idx2,idx3)
        B101(:) = B_grid(:,idx1p,idx2,idx3p)
        B110(:) = B_grid(:,idx1p,idx2p,idx3)
        B111(:) = B_grid(:,idx1p,idx2p,idx3p)

        delta1 = grid1(idx1p) - grid1(idx1)
        delta2 = grid2(idx2p) - grid2(idx2)
        delta3 = grid3(idx3p) - grid3(idx3)
        shift1 = (pos_in(1) - grid1(idx1))/delta1
        shift2 = (pos_in(2) - grid2(idx2))/delta2
        shift3 = (pos_in(3) - grid3(idx3))/delta3
        a1(:) = B100(:)-B000(:)
        a2(:) = B010(:)-B000(:)
        a3(:) = B001(:)-B000(:)
        a4(:) = B110(:)-B100(:)-B010(:)+B000(:)
        a5(:) = B101(:)-B100(:)-B001(:)+B000(:)
        a6(:) = B011(:)-B001(:)-B010(:)+B000(:)
        a7(:) = B111(:)-B110(:)-B101(:)-B011(:)+B001(:)+B010(:)+B100(:)-B000(:)

        B_out(:) = B000(:) + a1(:)*shift1 + a2(:)*shift2 + a3(:)*shift3 + a4(:)*shift1*shift2 &
               + a5(:)*shift1*shift3 + a6(:)*shift2*shift3 + a7(:)*shift1*shift2*shift3

      end subroutine B_interp_regular


      subroutine interp_single_regular(idx1, idx2, idx3, sz_1, sz_2, sz_3, grid1, grid2, grid3, &
                                        array_in, pos_in, value_out)
      !Trilinear interpolation of a single variable at a coordinate surrounded by grid cells
      !First transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]

        INTEGER :: idx1, idx2, idx3
        INTEGER :: sz_1, sz_2, sz_3
        REAL(num), DIMENSION(:), ALLOCATABLE :: grid1
        REAL(num), DIMENSION(:), ALLOCATABLE :: grid2
        REAL(num), DIMENSION(:), ALLOCATABLE :: grid3
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array_in
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: value_out

        INTEGER :: idx1p, idx2p, idx3p
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1, a2, a3, a4, a5, a6, a7
        REAL(num) :: arr000, arr001, arr010, arr011, arr100, arr101, arr110, arr111

        idx1 = find_index(grid1,sz_1,pos_in(1),idx1)
        idx2 = find_index(grid2,sz_2,pos_in(2),idx2)
        idx3 = find_index(grid3,sz_3,pos_in(3),idx3)
        idx1p = idx1+1
        idx2p = idx2+1
        idx3p = idx3+1
        arr000 = array_in(idx1,idx2,idx3)
        arr001 = array_in(idx1,idx2,idx3p)
        arr010 = array_in(idx1,idx2p,idx3)
        arr011 = array_in(idx1,idx2p,idx3p)
        arr100 = array_in(idx1p,idx2,idx3)
        arr101 = array_in(idx1p,idx2,idx3p)
        arr110 = array_in(idx1p,idx2p,idx3)
        arr111 = array_in(idx1p,idx2p,idx3p)

        delta1 = grid1(idx1p) - grid1(idx1)
        delta2 = grid2(idx2p) - grid2(idx2)
        delta3 = grid3(idx3p) - grid3(idx3)
        shift1 = (pos_in(1) - grid1(idx1))/delta1
        shift2 = (pos_in(2) - grid2(idx2))/delta2
        shift3 = (pos_in(3) - grid3(idx3))/delta3
        a1 = arr100-arr000
        a2 = arr010-arr000
        a3 = arr001-arr000
        a4 = arr110-arr100-arr010+arr000
        a5 = arr101-arr100-arr001+arr000
        a6 = arr011-arr001-arr010+arr000
        a7 = arr111-arr110-arr101-arr011+arr001+arr010+arr100-arr000

        value_out = arr000 + a1*shift1 + a2*shift2 + a3*shift3 + a4*shift1*shift2 &
               + a5*shift1*shift3 + a6*shift2*shift3 + a7*shift1*shift2*shift3

      end subroutine interp_single_regular


      subroutine B_interp_regular_separate(idx_in, pos_in, B_out)
      !Trilinear interpolation of B at a coordinate surrounded by grid cells
      !Each component of B has a separate array (e.g. staggered grid)
      !First transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)

        call interp_single_regular(idx_in(1), idx_in(2), idx_in(3), sz_11, sz_12, sz_13, grid1_1, &
                                        grid1_2, grid1_3, B_grid1, pos_in, B_out(1))
        call interp_single_regular(idx_in(4), idx_in(5), idx_in(6), sz_21, sz_22, sz_23, grid2_1, &
                                        grid2_2, grid2_3, B_grid2, pos_in, B_out(2))
        call interp_single_regular(idx_in(7), idx_in(8), idx_in(9), sz_31, sz_32, sz_33, grid3_1, &
                                        grid3_2, grid3_3, B_grid3, pos_in, B_out(3))

      end subroutine B_interp_regular_separate


      subroutine B_interp_normalized_regular(idx_in, pos_in, B_out)
      !Trilinear interpolation of B at a coordinate surrounded by grid cells
      !First normalize B at each vertex
      !First transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)

        INTEGER :: idx1, idx2, idx3
        INTEGER :: idx1p, idx2p, idx3p
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1(3), a2(3), a3(3), a4(3), a5(3), a6(3), a7(3)
        REAL(num) :: B000(3), B001(3), B010(3), B011(3), B100(3), B101(3), B110(3), B111(3)

        idx1 = idx_in(1)
        idx2 = idx_in(2)
        idx3 = idx_in(3)
        idx1 = find_index(grid1,sz_1,pos_in(1),idx1)
        idx2 = find_index(grid2,sz_2,pos_in(2),idx2)
        idx3 = find_index(grid3,sz_3,pos_in(3),idx3)
        idx_in(1) = idx1
        idx_in(2) = idx2
        idx_in(3) = idx3
        idx1p = idx1+1
        idx2p = idx2+1
        idx3p = idx3+1
        B000(:) = normalize_vector(B_grid(:,idx1,idx2,idx3))
        B001(:) = normalize_vector(B_grid(:,idx1,idx2,idx3p))
        B010(:) = normalize_vector(B_grid(:,idx1,idx2p,idx3))
        B011(:) = normalize_vector(B_grid(:,idx1,idx2p,idx3p))
        B100(:) = normalize_vector(B_grid(:,idx1p,idx2,idx3))
        B101(:) = normalize_vector(B_grid(:,idx1p,idx2,idx3p))
        B110(:) = normalize_vector(B_grid(:,idx1p,idx2p,idx3))
        B111(:) = normalize_vector(B_grid(:,idx1p,idx2p,idx3p))

        delta1 = grid1(idx1p) - grid1(idx1)
        delta2 = grid2(idx2p) - grid2(idx2)
        delta3 = grid3(idx3p) - grid3(idx3)
        shift1 = (pos_in(1) - grid1(idx1))/delta1
        shift2 = (pos_in(2) - grid2(idx2))/delta2
        shift3 = (pos_in(3) - grid3(idx3))/delta3
        a1(:) = B100(:)-B000(:)
        a2(:) = B010(:)-B000(:)
        a3(:) = B001(:)-B000(:)
        a4(:) = B110(:)-B100(:)-B010(:)+B000(:)
        a5(:) = B101(:)-B100(:)-B001(:)+B000(:)
        a6(:) = B011(:)-B001(:)-B010(:)+B000(:)
        a7(:) = B111(:)-B110(:)-B101(:)-B011(:)+B001(:)+B010(:)+B100(:)-B000(:)

        B_out(:) = B000(:) + a1(:)*shift1 + a2(:)*shift2 + a3(:)*shift3 + a4(:)*shift1*shift2 &
               + a5(:)*shift1*shift3 + a6(:)*shift2*shift3 + a7(:)*shift1*shift2*shift3

      end subroutine B_interp_normalized_regular


      subroutine B_interp_irregular(idx_in, pos_in, B_out)
      !Trilinear interpolation of B at a coordinate surrounded by grid cells
      !First transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)

        INTEGER :: idx_b
        INTEGER :: idx1m, idx2m, idx3m, idx1, idx2, idx3, idx1p, idx2p, idx3p
        REAL(num) :: coord1l, coord1r, coord2l, coord2r, coord3l, coord3r
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1(3), a2(3), a3(3), a4(3), a5(3), a6(3), a7(3)
        REAL(num) :: B000(3), B001(3), B010(3), B011(3), B100(3), B101(3), B110(3), B111(3)

        idx_b = idx_in(1)
        idx_b = find_index_irregular(pos_in,idx_b)
        idx_in(1) = idx_b
        idx1m = FLOOR((pos_in(1)-grid1_ir(1,idx_b))/(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                     *REAL(sz_1-1,num))
        if (idx1m .ge. (sz_1-1)) then
          idx1m = sz_1 - 2
        end if
        idx2m = FLOOR((pos_in(2)-grid2_ir(1,idx_b))/(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                     *REAL(sz_2-1,num))
        if (idx2m .ge. (sz_2-1)) then
          idx2m = sz_2 - 2
        end if
        idx3m = FLOOR((pos_in(3)-grid3_ir(1,idx_b))/(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                     *REAL(sz_3-1,num))
        if (idx3m .ge. (sz_3-1)) then
          idx3m = sz_3 - 2
        end if
        idx1 = idx1m+1
        idx2 = idx2m+1
        idx3 = idx3m+1
        idx1p = idx1m+2
        idx2p = idx2m+2
        idx3p = idx3m+2
        B000(:) = B_grid_ir(:,idx1,idx2,idx3,idx_b)
        B001(:) = B_grid_ir(:,idx1,idx2,idx3p,idx_b)
        B010(:) = B_grid_ir(:,idx1,idx2p,idx3,idx_b)
        B011(:) = B_grid_ir(:,idx1,idx2p,idx3p,idx_b)
        B100(:) = B_grid_ir(:,idx1p,idx2,idx3,idx_b)
        B101(:) = B_grid_ir(:,idx1p,idx2,idx3p,idx_b)
        B110(:) = B_grid_ir(:,idx1p,idx2p,idx3,idx_b)
        B111(:) = B_grid_ir(:,idx1p,idx2p,idx3p,idx_b)

        coord1l = REAL(idx1m,num)/REAL(sz_1-1,num)*(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                  +grid1_ir(1,idx_b)
        coord1r = REAL(idx1,num)/REAL(sz_1-1,num)*(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                  +grid1_ir(1,idx_b)
        coord2l = REAL(idx2m,num)/REAL(sz_2-1,num)*(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                  +grid2_ir(1,idx_b)
        coord2r = REAL(idx2,num)/REAL(sz_2-1,num)*(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                  +grid2_ir(1,idx_b)
        coord3l = REAL(idx3m,num)/REAL(sz_3-1,num)*(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                  +grid3_ir(1,idx_b)
        coord3r = REAL(idx3,num)/REAL(sz_3-1,num)*(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                  +grid3_ir(1,idx_b)
        delta1 = coord1r - coord1l
        delta2 = coord2r - coord2l
        delta3 = coord3r - coord3l
        shift1 = (pos_in(1) - coord1l)/delta1
        shift2 = (pos_in(2) - coord2l)/delta2
        shift3 = (pos_in(3) - coord3l)/delta3
        a1(:) = B100(:)-B000(:)
        a2(:) = B010(:)-B000(:)
        a3(:) = B001(:)-B000(:)
        a4(:) = B110(:)-B100(:)-B010(:)+B000(:)
        a5(:) = B101(:)-B100(:)-B001(:)+B000(:)
        a6(:) = B011(:)-B001(:)-B010(:)+B000(:)
        a7(:) = B111(:)-B110(:)-B101(:)-B011(:)+B001(:)+B010(:)+B100(:)-B000(:)

        B_out(:) = B000(:) + a1(:)*shift1 + a2(:)*shift2 + a3(:)*shift3 + a4(:)*shift1*shift2 &
               + a5(:)*shift1*shift3 + a6(:)*shift2*shift3 + a7(:)*shift1*shift2*shift3

      end subroutine B_interp_irregular


      subroutine B_interp_normalized_irregular(idx_in, pos_in, B_out)
      !Trilinear interpolation of B at a coordinate surrounded by grid cells
      !First normalize B at each vertex
      !First transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)

        INTEGER :: idx_b
        INTEGER :: idx1m, idx2m, idx3m, idx1, idx2, idx3, idx1p, idx2p, idx3p
        REAL(num) :: coord1l, coord1r, coord2l, coord2r, coord3l, coord3r
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1(3), a2(3), a3(3), a4(3), a5(3), a6(3), a7(3)
        REAL(num) :: B000(3), B001(3), B010(3), B011(3), B100(3), B101(3), B110(3), B111(3)

        idx_b = idx_in(1)
        idx_b = find_index_irregular(pos_in,idx_b)
        idx_in(1) = idx_b
        idx1m = FLOOR((pos_in(1)-grid1_ir(1,idx_b))/(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                     *REAL(sz_1-1,num))
        if (idx1m .ge. (sz_1-1)) then
          idx1m = sz_1 - 2
        end if
        idx2m = FLOOR((pos_in(2)-grid2_ir(1,idx_b))/(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                     *REAL(sz_2-1,num))
        if (idx2m .ge. (sz_2-1)) then
          idx2m = sz_2 - 2
        end if
        idx3m = FLOOR((pos_in(3)-grid3_ir(1,idx_b))/(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                     *REAL(sz_3-1,num))
        if (idx3m .ge. (sz_3-1)) then
          idx3m = sz_3 - 2
        end if
        idx1 = idx1m+1
        idx2 = idx2m+1
        idx3 = idx3m+1
        idx1p = idx1m+2
        idx2p = idx2m+2
        idx3p = idx3m+2
        B000(:) = normalize_vector(B_grid_ir(:,idx1,idx2,idx3,idx_b))
        B001(:) = normalize_vector(B_grid_ir(:,idx1,idx2,idx3p,idx_b))
        B010(:) = normalize_vector(B_grid_ir(:,idx1,idx2p,idx3,idx_b))
        B011(:) = normalize_vector(B_grid_ir(:,idx1,idx2p,idx3p,idx_b))
        B100(:) = normalize_vector(B_grid_ir(:,idx1p,idx2,idx3,idx_b))
        B101(:) = normalize_vector(B_grid_ir(:,idx1p,idx2,idx3p,idx_b))
        B110(:) = normalize_vector(B_grid_ir(:,idx1p,idx2p,idx3,idx_b))
        B111(:) = normalize_vector(B_grid_ir(:,idx1p,idx2p,idx3p,idx_b))

        coord1l = REAL(idx1m,num)/REAL(sz_1-1,num)*(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                  +grid1_ir(1,idx_b)
        coord1r = REAL(idx1,num)/REAL(sz_1-1,num)*(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                  +grid1_ir(1,idx_b)
        coord2l = REAL(idx2m,num)/REAL(sz_2-1,num)*(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                  +grid2_ir(1,idx_b)
        coord2r = REAL(idx2,num)/REAL(sz_2-1,num)*(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                  +grid2_ir(1,idx_b)
        coord3l = REAL(idx3m,num)/REAL(sz_3-1,num)*(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                  +grid3_ir(1,idx_b)
        coord3r = REAL(idx3,num)/REAL(sz_3-1,num)*(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                  +grid3_ir(1,idx_b)
        delta1 = coord1r - coord1l
        delta2 = coord2r - coord2l
        delta3 = coord3r - coord3l
        shift1 = (pos_in(1) - coord1l)/delta1
        shift2 = (pos_in(2) - coord2l)/delta2
        shift3 = (pos_in(3) - coord3l)/delta3
        a1(:) = B100(:)-B000(:)
        a2(:) = B010(:)-B000(:)
        a3(:) = B001(:)-B000(:)
        a4(:) = B110(:)-B100(:)-B010(:)+B000(:)
        a5(:) = B101(:)-B100(:)-B001(:)+B000(:)
        a6(:) = B011(:)-B001(:)-B010(:)+B000(:)
        a7(:) = B111(:)-B110(:)-B101(:)-B011(:)+B001(:)+B010(:)+B100(:)-B000(:)

        B_out(:) = B000(:) + a1(:)*shift1 + a2(:)*shift2 + a3(:)*shift3 + a4(:)*shift1*shift2 &
               + a5(:)*shift1*shift3 + a6(:)*shift2*shift3 + a7(:)*shift1*shift2*shift3

      end subroutine B_interp_normalized_irregular


      subroutine B_gradB_interp_regular(idx_in, pos_in, B_out, gradB_out)
      !Trilinear interpolation of B at a coordinate surrounded by grid cells
      !Transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]
      !Also returns gradient, gradB_out(i,j) is d B_i / d x_j
      !Note: Cartesian only, needs to be divided by R etc afterwards for spherical (effectively just partial derivatives)

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)
        REAL(num), INTENT(OUT) :: gradB_out(3,3)

        INTEGER :: idx1, idx2, idx3
        INTEGER :: idx1p, idx2p, idx3p
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1(3), a2(3), a3(3), a4(3), a5(3), a6(3), a7(3)
        REAL(num) :: B000(3), B001(3), B010(3), B011(3), B100(3), B101(3), B110(3), B111(3)

        idx1 = idx_in(1)
        idx2 = idx_in(2)
        idx3 = idx_in(3)
        idx1 = find_index(grid1,sz_1,pos_in(1),idx1)
        idx2 = find_index(grid2,sz_2,pos_in(2),idx2)
        idx3 = find_index(grid3,sz_3,pos_in(3),idx3)
        idx_in(1) = idx1
        idx_in(2) = idx2
        idx_in(3) = idx3
        idx1p = idx1+1
        idx2p = idx2+1
        idx3p = idx3+1
        B000(:) = B_grid(:,idx1,idx2,idx3)
        B001(:) = B_grid(:,idx1,idx2,idx3p)
        B010(:) = B_grid(:,idx1,idx2p,idx3)
        B011(:) = B_grid(:,idx1,idx2p,idx3p)
        B100(:) = B_grid(:,idx1p,idx2,idx3)
        B101(:) = B_grid(:,idx1p,idx2,idx3p)
        B110(:) = B_grid(:,idx1p,idx2p,idx3)
        B111(:) = B_grid(:,idx1p,idx2p,idx3p)

        delta1 = grid1(idx1p) - grid1(idx1)
        delta2 = grid2(idx2p) - grid2(idx2)
        delta3 = grid3(idx3p) - grid3(idx3)
        shift1 = (pos_in(1) - grid1(idx1))/delta1
        shift2 = (pos_in(2) - grid2(idx2))/delta2
        shift3 = (pos_in(3) - grid3(idx3))/delta3
        a1(:) = B100(:)-B000(:)
        a2(:) = B010(:)-B000(:)
        a3(:) = B001(:)-B000(:)
        a4(:) = B110(:)-B100(:)-B010(:)+B000(:)
        a5(:) = B101(:)-B100(:)-B001(:)+B000(:)
        a6(:) = B011(:)-B001(:)-B010(:)+B000(:)
        a7(:) = B111(:)-B110(:)-B101(:)-B011(:)+B001(:)+B010(:)+B100(:)-B000(:)

        B_out(:) = B000(:) + a1(:)*shift1 + a2(:)*shift2 + a3(:)*shift3 + a4(:)*shift1*shift2 &
               + a5(:)*shift1*shift3 + a6(:)*shift2*shift3 + a7(:)*shift1*shift2*shift3

        gradB_out(:,1) = (a1(:) + a4(:)*shift2 + a5(:)*shift3 + a7(:)*shift2*shift3) / delta1
        gradB_out(:,2) = (a2(:) + a4(:)*shift1 + a6(:)*shift3 + a7(:)*shift1*shift3) / delta2
        gradB_out(:,3) = (a3(:) + a5(:)*shift1 + a6(:)*shift2 + a7(:)*shift1*shift2) / delta3

      end subroutine B_gradB_interp_regular


      subroutine grad_interp_single_regular(idx1, idx2, idx3, sz_1, sz_2, sz_3, grid1, &
                                            grid2, grid3, array_in, pos_in, value_out, grad_out)
      !Trilinear interpolation of B at a coordinate surrounded by grid cells
      !Transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]
      !Also returns gradient, grad_out(j) is d array_in / d x_j
      !Note: Cartesian only, needs to be divided by R etc afterwards for spherical (effectively just partial derivatives)

        INTEGER :: idx1, idx2, idx3
        INTEGER :: sz_1, sz_2, sz_3
        REAL(num), DIMENSION(:), ALLOCATABLE :: grid1
        REAL(num), DIMENSION(:), ALLOCATABLE :: grid2
        REAL(num), DIMENSION(:), ALLOCATABLE :: grid3
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array_in
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: value_out
        REAL(num), INTENT(OUT) :: grad_out(3)

        INTEGER :: idx1p, idx2p, idx3p
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1, a2, a3, a4, a5, a6, a7
        REAL(num) :: arr000, arr001, arr010, arr011, arr100, arr101, arr110, arr111

        idx1 = find_index(grid1,sz_1,pos_in(1),idx1)
        idx2 = find_index(grid2,sz_2,pos_in(2),idx2)
        idx3 = find_index(grid3,sz_3,pos_in(3),idx3)
        idx1p = idx1+1
        idx2p = idx2+1
        idx3p = idx3+1
        arr000 = array_in(idx1,idx2,idx3)
        arr001 = array_in(idx1,idx2,idx3p)
        arr010 = array_in(idx1,idx2p,idx3)
        arr011 = array_in(idx1,idx2p,idx3p)
        arr100 = array_in(idx1p,idx2,idx3)
        arr101 = array_in(idx1p,idx2,idx3p)
        arr110 = array_in(idx1p,idx2p,idx3)
        arr111 = array_in(idx1p,idx2p,idx3p)

        delta1 = grid1(idx1p) - grid1(idx1)
        delta2 = grid2(idx2p) - grid2(idx2)
        delta3 = grid3(idx3p) - grid3(idx3)
        shift1 = (pos_in(1) - grid1(idx1))/delta1
        shift2 = (pos_in(2) - grid2(idx2))/delta2
        shift3 = (pos_in(3) - grid3(idx3))/delta3
        a1 = arr100-arr000
        a2 = arr010-arr000
        a3 = arr001-arr000
        a4 = arr110-arr100-arr010+arr000
        a5 = arr101-arr100-arr001+arr000
        a6 = arr011-arr001-arr010+arr000
        a7 = arr111-arr110-arr101-arr011+arr001+arr010+arr100-arr000

        value_out = arr000 + a1*shift1 + a2*shift2 + a3*shift3 + a4*shift1*shift2 &
               + a5*shift1*shift3 + a6*shift2*shift3 + a7*shift1*shift2*shift3

        grad_out(1) = (a1 + a4*shift2 + a5*shift3 + a7*shift2*shift3) / delta1
        grad_out(2) = (a2 + a4*shift1 + a6*shift3 + a7*shift1*shift3) / delta2
        grad_out(3) = (a3 + a5*shift1 + a6*shift2 + a7*shift1*shift2) / delta3

      end subroutine grad_interp_single_regular


      subroutine B_gradB_interp_regular_separate(idx_in, pos_in, B_out, gradB_out)
      !Trilinear interpolation of B at a coordinate surrounded by grid cells
      !Each component of B has a separate array (e.g. staggered grid)
      !First transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)
        REAL(num), INTENT(OUT) :: gradB_out(3,3)

        REAL(num) :: grad_temp(3)

        call grad_interp_single_regular(idx_in(1), idx_in(2), idx_in(3), sz_11, sz_12, sz_13, &
                              grid1_1, grid1_2, grid1_3, B_grid1, pos_in, B_out(1), grad_temp)
        gradB_out(1,:) = grad_temp(:)
        call grad_interp_single_regular(idx_in(4), idx_in(5), idx_in(6), sz_21, sz_22, sz_23, &
                              grid2_1, grid2_2, grid2_3, B_grid2, pos_in, B_out(2), grad_temp)
        gradB_out(2,:) = grad_temp(:)
        call grad_interp_single_regular(idx_in(7), idx_in(8), idx_in(9), sz_31, sz_32, sz_33, &
                              grid3_1, grid3_2, grid3_3, B_grid3, pos_in, B_out(3), grad_temp)
        gradB_out(3,:) = grad_temp(:)

      end subroutine B_gradB_interp_regular_separate


      subroutine B_gradB_interp_normalized_regular(idx_in, pos_in, B_out, gradB_out)
      !Trilinear interpolation of norm(B) at a coordinate surrounded by grid cells
      !First normalize B at each vertex
      !Transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]
      !Also returns gradient, gradB_out(i,j) is d B_i / d x_j
      !Note: Cartesian only, needs to be divided by R etc afterwards for spherical (effectively just partial derivatives)

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)
        REAL(num), INTENT(OUT) :: gradB_out(3,3)

        INTEGER :: idx1, idx2, idx3
        INTEGER :: idx1p, idx2p, idx3p
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1(3), a2(3), a3(3), a4(3), a5(3), a6(3), a7(3)
        REAL(num) :: B000(3), B001(3), B010(3), B011(3), B100(3), B101(3), B110(3), B111(3)

        idx1 = idx_in(1)
        idx2 = idx_in(2)
        idx3 = idx_in(3)
        idx1 = find_index(grid1,sz_1,pos_in(1),idx1)
        idx2 = find_index(grid2,sz_2,pos_in(2),idx2)
        idx3 = find_index(grid3,sz_3,pos_in(3),idx3)
        idx_in(1) = idx1
        idx_in(2) = idx2
        idx_in(3) = idx3
        idx1p = idx1+1
        idx2p = idx2+1
        idx3p = idx3+1
        B000(:) = normalize_vector(B_grid(:,idx1,idx2,idx3))
        B001(:) = normalize_vector(B_grid(:,idx1,idx2,idx3p))
        B010(:) = normalize_vector(B_grid(:,idx1,idx2p,idx3))
        B011(:) = normalize_vector(B_grid(:,idx1,idx2p,idx3p))
        B100(:) = normalize_vector(B_grid(:,idx1p,idx2,idx3))
        B101(:) = normalize_vector(B_grid(:,idx1p,idx2,idx3p))
        B110(:) = normalize_vector(B_grid(:,idx1p,idx2p,idx3))
        B111(:) = normalize_vector(B_grid(:,idx1p,idx2p,idx3p))

        delta1 = grid1(idx1p) - grid1(idx1)
        delta2 = grid2(idx2p) - grid2(idx2)
        delta3 = grid3(idx3p) - grid3(idx3)
        shift1 = (pos_in(1) - grid1(idx1))/delta1
        shift2 = (pos_in(2) - grid2(idx2))/delta2
        shift3 = (pos_in(3) - grid3(idx3))/delta3
        a1(:) = B100(:)-B000(:)
        a2(:) = B010(:)-B000(:)
        a3(:) = B001(:)-B000(:)
        a4(:) = B110(:)-B100(:)-B010(:)+B000(:)
        a5(:) = B101(:)-B100(:)-B001(:)+B000(:)
        a6(:) = B011(:)-B001(:)-B010(:)+B000(:)
        a7(:) = B111(:)-B110(:)-B101(:)-B011(:)+B001(:)+B010(:)+B100(:)-B000(:)

        B_out(:) = B000(:) + a1(:)*shift1 + a2(:)*shift2 + a3(:)*shift3 + a4(:)*shift1*shift2 &
               + a5(:)*shift1*shift3 + a6(:)*shift2*shift3 + a7(:)*shift1*shift2*shift3

        gradB_out(:,1) = (a1(:) + a4(:)*shift2 + a5(:)*shift3 + a7(:)*shift2*shift3) / delta1
        gradB_out(:,2) = (a2(:) + a4(:)*shift1 + a6(:)*shift3 + a7(:)*shift1*shift3) / delta2
        gradB_out(:,3) = (a3(:) + a5(:)*shift1 + a6(:)*shift2 + a7(:)*shift1*shift2) / delta3

      end subroutine B_gradB_interp_normalized_regular


      subroutine B_gradB_interp_irregular(idx_in, pos_in, B_out, gradB_out)
      !Trilinear interpolation of norm(B) at a coordinate surrounded by grid cells
      !First normalize B at each vertex
      !Transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]
      !Also returns gradient, gradB_out(i,j) is d B_i / d x_j
      !Note: Cartesian only, needs to be divided by R etc afterwards for spherical (effectively just partial derivatives)

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)
        REAL(num), INTENT(OUT) :: gradB_out(3,3)

        INTEGER :: idx_b
        INTEGER :: idx1m, idx2m, idx3m, idx1, idx2, idx3, idx1p, idx2p, idx3p
        REAL(num) :: coord1l, coord1r, coord2l, coord2r, coord3l, coord3r
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1(3), a2(3), a3(3), a4(3), a5(3), a6(3), a7(3)
        REAL(num) :: B000(3), B001(3), B010(3), B011(3), B100(3), B101(3), B110(3), B111(3)

        idx_b = idx_in(1)
        idx_b = find_index_irregular(pos_in,idx_b)
        idx_in(1) = idx_b
        idx1m = FLOOR((pos_in(1)-grid1_ir(1,idx_b))/(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                     *REAL(sz_1-1,num))
        if (idx1m .ge. (sz_1-1)) then
          idx1m = sz_1 - 2
        end if
        idx2m = FLOOR((pos_in(2)-grid2_ir(1,idx_b))/(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                     *REAL(sz_2-1,num))
        if (idx2m .ge. (sz_2-1)) then
          idx2m = sz_2 - 2
        end if
        idx3m = FLOOR((pos_in(3)-grid3_ir(1,idx_b))/(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                     *REAL(sz_3-1,num))
        if (idx3m .ge. (sz_3-1)) then
          idx3m = sz_3 - 2
        end if
        idx1 = idx1m+1
        idx2 = idx2m+1
        idx3 = idx3m+1
        idx1p = idx1m+2
        idx2p = idx2m+2
        idx3p = idx3m+2
        B000(:) = normalize_vector(B_grid_ir(:,idx1,idx2,idx3,idx_b))
        B001(:) = normalize_vector(B_grid_ir(:,idx1,idx2,idx3p,idx_b))
        B010(:) = normalize_vector(B_grid_ir(:,idx1,idx2p,idx3,idx_b))
        B011(:) = normalize_vector(B_grid_ir(:,idx1,idx2p,idx3p,idx_b))
        B100(:) = normalize_vector(B_grid_ir(:,idx1p,idx2,idx3,idx_b))
        B101(:) = normalize_vector(B_grid_ir(:,idx1p,idx2,idx3p,idx_b))
        B110(:) = normalize_vector(B_grid_ir(:,idx1p,idx2p,idx3,idx_b))
        B111(:) = normalize_vector(B_grid_ir(:,idx1p,idx2p,idx3p,idx_b))

        coord1l = REAL(idx1m,num)/REAL(sz_1-1,num)*(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                  +grid1_ir(1,idx_b)
        coord1r = REAL(idx1,num)/REAL(sz_1-1,num)*(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                  +grid1_ir(1,idx_b)
        coord2l = REAL(idx2m,num)/REAL(sz_2-1,num)*(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                  +grid2_ir(1,idx_b)
        coord2r = REAL(idx2,num)/REAL(sz_2-1,num)*(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                  +grid2_ir(1,idx_b)
        coord3l = REAL(idx3m,num)/REAL(sz_3-1,num)*(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                  +grid3_ir(1,idx_b)
        coord3r = REAL(idx3,num)/REAL(sz_3-1,num)*(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                  +grid3_ir(1,idx_b)
        delta1 = coord1r - coord1l
        delta2 = coord2r - coord2l
        delta3 = coord3r - coord3l
        shift1 = (pos_in(1) - coord1l)/delta1
        shift2 = (pos_in(2) - coord2l)/delta2
        shift3 = (pos_in(3) - coord3l)/delta3
        a1(:) = B100(:)-B000(:)
        a2(:) = B010(:)-B000(:)
        a3(:) = B001(:)-B000(:)
        a4(:) = B110(:)-B100(:)-B010(:)+B000(:)
        a5(:) = B101(:)-B100(:)-B001(:)+B000(:)
        a6(:) = B011(:)-B001(:)-B010(:)+B000(:)
        a7(:) = B111(:)-B110(:)-B101(:)-B011(:)+B001(:)+B010(:)+B100(:)-B000(:)

        B_out(:) = B000(:) + a1(:)*shift1 + a2(:)*shift2 + a3(:)*shift3 + a4(:)*shift1*shift2 &
               + a5(:)*shift1*shift3 + a6(:)*shift2*shift3 + a7(:)*shift1*shift2*shift3

        gradB_out(:,1) = (a1(:) + a4(:)*shift2 + a5(:)*shift3 + a7(:)*shift2*shift3) / delta1
        gradB_out(:,2) = (a2(:) + a4(:)*shift1 + a6(:)*shift3 + a7(:)*shift1*shift3) / delta2
        gradB_out(:,3) = (a3(:) + a5(:)*shift1 + a6(:)*shift2 + a7(:)*shift1*shift2) / delta3

      end subroutine B_gradB_interp_irregular


      subroutine B_gradB_interp_normalized_irregular(idx_in, pos_in, B_out, gradB_out)
      !Trilinear interpolation of norm(B) at a coordinate surrounded by grid cells
      !First normalize B at each vertex
      !Transform to a cell of size [0,1]
      !delta1,2,3 are cell sizes
      !shift1,2,3 are the coordinates transformed to [0,1]
      !Also returns gradient, gradB_out(i,j) is d B_i / d x_j
      !Note: Cartesian only, needs to be divided by R etc afterwards for spherical (effectively just partial derivatives)

        INTEGER :: idx_in(9)
        REAL(num), INTENT(IN) :: pos_in(3)
        REAL(num), INTENT(OUT) :: B_out(3)
        REAL(num), INTENT(OUT) :: gradB_out(3,3)

        INTEGER :: idx_b
        INTEGER :: idx1m, idx2m, idx3m, idx1, idx2, idx3, idx1p, idx2p, idx3p
        REAL(num) :: coord1l, coord1r, coord2l, coord2r, coord3l, coord3r
        REAL(num) :: shift1, shift2, shift3, delta1, delta2, delta3
        REAL(num) :: a1(3), a2(3), a3(3), a4(3), a5(3), a6(3), a7(3)
        REAL(num) :: B000(3), B001(3), B010(3), B011(3), B100(3), B101(3), B110(3), B111(3)

        idx_b = idx_in(1)
        idx_b = find_index_irregular(pos_in,idx_b)
        idx_in(1) = idx_b
        idx1m = FLOOR((pos_in(1)-grid1_ir(1,idx_b))/(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                     *REAL(sz_1-1,num))
        if (idx1m .ge. (sz_1-1)) then
          idx1m = sz_1 - 2
        end if
        idx2m = FLOOR((pos_in(2)-grid2_ir(1,idx_b))/(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                     *REAL(sz_2-1,num))
        if (idx2m .ge. (sz_2-1)) then
          idx2m = sz_2 - 2
        end if
        idx3m = FLOOR((pos_in(3)-grid3_ir(1,idx_b))/(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                     *REAL(sz_3-1,num))
        if (idx3m .ge. (sz_3-1)) then
          idx3m = sz_3 - 2
        end if
        idx1 = idx1m+1
        idx2 = idx2m+1
        idx3 = idx3m+1
        idx1p = idx1m+2
        idx2p = idx2m+2
        idx3p = idx3m+2
        B000(:) = normalize_vector(B_grid_ir(:,idx1,idx2,idx3,idx_b))
        B001(:) = normalize_vector(B_grid_ir(:,idx1,idx2,idx3p,idx_b))
        B010(:) = normalize_vector(B_grid_ir(:,idx1,idx2p,idx3,idx_b))
        B011(:) = normalize_vector(B_grid_ir(:,idx1,idx2p,idx3p,idx_b))
        B100(:) = normalize_vector(B_grid_ir(:,idx1p,idx2,idx3,idx_b))
        B101(:) = normalize_vector(B_grid_ir(:,idx1p,idx2,idx3p,idx_b))
        B110(:) = normalize_vector(B_grid_ir(:,idx1p,idx2p,idx3,idx_b))
        B111(:) = normalize_vector(B_grid_ir(:,idx1p,idx2p,idx3p,idx_b))

        coord1l = REAL(idx1m,num)/REAL(sz_1-1,num)*(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                  +grid1_ir(1,idx_b)
        coord1r = REAL(idx1,num)/REAL(sz_1-1,num)*(grid1_ir(2,idx_b)-grid1_ir(1,idx_b)) &
                  +grid1_ir(1,idx_b)
        coord2l = REAL(idx2m,num)/REAL(sz_2-1,num)*(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                  +grid2_ir(1,idx_b)
        coord2r = REAL(idx2,num)/REAL(sz_2-1,num)*(grid2_ir(2,idx_b)-grid2_ir(1,idx_b)) &
                  +grid2_ir(1,idx_b)
        coord3l = REAL(idx3m,num)/REAL(sz_3-1,num)*(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                  +grid3_ir(1,idx_b)
        coord3r = REAL(idx3,num)/REAL(sz_3-1,num)*(grid3_ir(2,idx_b)-grid3_ir(1,idx_b)) &
                  +grid3_ir(1,idx_b)
        delta1 = coord1r - coord1l
        delta2 = coord2r - coord2l
        delta3 = coord3r - coord3l
        shift1 = (pos_in(1) - coord1l)/delta1
        shift2 = (pos_in(2) - coord2l)/delta2
        shift3 = (pos_in(3) - coord3l)/delta3
        a1(:) = B100(:)-B000(:)
        a2(:) = B010(:)-B000(:)
        a3(:) = B001(:)-B000(:)
        a4(:) = B110(:)-B100(:)-B010(:)+B000(:)
        a5(:) = B101(:)-B100(:)-B001(:)+B000(:)
        a6(:) = B011(:)-B001(:)-B010(:)+B000(:)
        a7(:) = B111(:)-B110(:)-B101(:)-B011(:)+B001(:)+B010(:)+B100(:)-B000(:)

        B_out(:) = B000(:) + a1(:)*shift1 + a2(:)*shift2 + a3(:)*shift3 + a4(:)*shift1*shift2 &
               + a5(:)*shift1*shift3 + a6(:)*shift2*shift3 + a7(:)*shift1*shift2*shift3

        gradB_out(:,1) = (a1(:) + a4(:)*shift2 + a5(:)*shift3 + a7(:)*shift2*shift3) / delta1
        gradB_out(:,2) = (a2(:) + a4(:)*shift1 + a6(:)*shift3 + a7(:)*shift1*shift3) / delta2
        gradB_out(:,3) = (a3(:) + a5(:)*shift1 + a6(:)*shift2 + a7(:)*shift1*shift2) / delta3

      end subroutine B_gradB_interp_normalized_irregular


      subroutine step_cartesian(idx_in, pos_in, pos_out, u_vec, v_vec, dl, &
                                keep_running, check_position, B_interp, B_gradB_interp)

        INTEGER :: idx_in(9)
        REAL(num) :: pos_in(3),pos_out(3),u_vec(3),v_vec(3),dl
        LOGICAL :: keep_running
        procedure(check_position_iface) :: check_position
        procedure(B_interp_iface) :: B_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp

        REAL(num) :: mod_B, pos_0out(3), k_1(3), B_curr(3), dl_norm

        call check_position(pos_in,pos_0out,keep_running)
        call B_interp(idx_in, pos_0out, B_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_1(:) = B_curr(:)*dl_norm

        pos_out(:) = pos_in(:) + k_1(:)
        call check_position(pos_out,pos_0out,keep_running)

      end subroutine step_cartesian


      subroutine step_cartesian_RK4(idx_in, pos_in, pos_out, dl, keep_running, &
                                check_position, B_interp, B_gradB_interp) !TODO - fix and use this

        INTEGER :: idx_in(9)
        REAL(num) :: pos_in(3),pos_out(3),dl
        LOGICAL :: keep_running
        procedure(check_position_iface) :: check_position
        procedure(B_interp_iface) :: B_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp

        REAL(num) :: mod_B, pos_1in(3), pos_2in(3), pos_3in(3), pos_0out(3), pos_1out(3), dl_norm
        REAL(num) :: pos_2out(3), pos_3out(3), pos_4out(3), k_1(3), k_2(3), k_3(3), k_4(3)
        REAL(num) :: B_curr(3)

        call check_position(pos_in,pos_0out,keep_running)
        call B_interp(idx_in, pos_0out, B_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_1(:) = B_curr(:)*dl_norm

        pos_1in(:) = pos_in(:) + 0.5_num*k_1(:) 
        call check_position(pos_1in,pos_1out,keep_running)
        call B_interp(idx_in, pos_1out, B_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_2(:) = B_curr(:)*dl_norm

        pos_2in(:) = pos_in(:) + 0.5_num*k_2(:)
        call check_position(pos_2in,pos_2out,keep_running)
        call B_interp(idx_in, pos_2out, B_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_3(:) = B_curr(:)*dl_norm

        pos_3in(:) = pos_in(:) + k_3(:)
        call check_position(pos_3in,pos_3out,keep_running)
        call B_interp(idx_in, pos_3out, B_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_4(:) = B_curr(:)*dl_norm

        pos_out(:) = pos_in(:) + (k_1(:) + 2.0_num*k_2(:) + 2.0_num*k_3(:) + k_4(:))/6.0_num
        call check_position(pos_out,pos_4out,keep_running)

      end subroutine step_cartesian_RK4


      subroutine step_cartesianQ(idx_in, pos_in, pos_out, u_vec, v_vec, dl, &
                                 keep_running, check_position, B_interp, B_gradB_interp)

        INTEGER :: idx_in(9)
        REAL(num) :: pos_in(3),pos_out(3),u_vec(3),v_vec(3),dl
        LOGICAL :: keep_running
        procedure(check_position_iface) :: check_position
        procedure(B_interp_iface) :: B_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp

        REAL(num) :: mod_B, pos_0out(3), k_1(3), B_curr(3), dl_norm
        REAL(num) :: k_1u(3), k_1v(3), gradB_curr(3,3)

        call check_position(pos_in,pos_0out,keep_running)
        call B_gradB_interp(idx_in, pos_0out, B_curr, gradB_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_1(:) = B_curr(:)*dl_norm
        k_1u(1) = (u_vec(1)*gradB_curr(1,1)+u_vec(2)*gradB_curr(1,2)+u_vec(3)*gradB_curr(1,3))&
                  *dl_norm
        k_1u(2) = (u_vec(1)*gradB_curr(2,1)+u_vec(2)*gradB_curr(2,2)+u_vec(3)*gradB_curr(2,3))&
                  *dl_norm
        k_1u(3) = (u_vec(1)*gradB_curr(3,1)+u_vec(2)*gradB_curr(3,2)+u_vec(3)*gradB_curr(3,3))&
                  *dl_norm
        k_1v(1) = (v_vec(1)*gradB_curr(1,1)+v_vec(2)*gradB_curr(1,2)+v_vec(3)*gradB_curr(1,3))&
                  *dl_norm
        k_1v(2) = (v_vec(1)*gradB_curr(2,1)+v_vec(2)*gradB_curr(2,2)+v_vec(3)*gradB_curr(2,3))&
                  *dl_norm
        k_1v(3) = (v_vec(1)*gradB_curr(3,1)+v_vec(2)*gradB_curr(3,2)+v_vec(3)*gradB_curr(3,3))&
                  *dl_norm

        pos_out(:) = pos_in(:) + k_1(:)
        call check_position(pos_out,pos_0out,keep_running)
        IF (keep_running) THEN
          u_vec(:) = u_vec(:) + k_1u(:)
          v_vec(:) = v_vec(:) + k_1v(:)
        END IF

      end subroutine step_cartesianQ


      subroutine step_cartesianQ_RK4(idx_in, pos_in, pos_out, u_vec, v_vec, dl, &
                                 keep_running, check_position, B_interp, B_gradB_interp) !TODO - fix and use this

        INTEGER :: idx_in(9)
        REAL(num) :: pos_in(3),pos_out(3),u_vec(3),v_vec(3),dl
        LOGICAL :: keep_running
        procedure(check_position_iface) :: check_position
        procedure(B_interp_iface) :: B_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp

        REAL(num) :: mod_B, dl_norm
        REAL(num) :: pos_1in(3), pos_2in(3), pos_3in(3), pos_0out(3), pos_1out(3), pos_2out(3)
        REAL(num) :: pos_3out(3), pos_4out(3), k_1(3), k_2(3), k_3(3), k_4(3), B_curr(3)
        REAL(num) :: u_1(3), u_2(3), u_3(3), v_1(3), v_2(3), v_3(3)
        REAL(num) :: k_1u(3), k_2u(3), k_3u(3), k_4u(3), k_1v(3), k_2v(3), k_3v(3), k_4v(3)
        REAL(num) :: gradB_curr(3,3)

        call check_position(pos_in,pos_0out,keep_running)
        call B_gradB_interp(idx_in, pos_0out, B_curr, gradB_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_1(:) = B_curr(:)*dl_norm
        k_1u(1) = (u_vec(1)*gradB_curr(1,1)+u_vec(2)*gradB_curr(1,2)+u_vec(3)*gradB_curr(1,3))&
                  *dl_norm
        k_1u(2) = (u_vec(1)*gradB_curr(2,1)+u_vec(2)*gradB_curr(2,2)+u_vec(3)*gradB_curr(2,3))&
                  *dl_norm
        k_1u(3) = (u_vec(1)*gradB_curr(3,1)+u_vec(2)*gradB_curr(3,2)+u_vec(3)*gradB_curr(3,3))&
                  *dl_norm
        k_1v(1) = (v_vec(1)*gradB_curr(1,1)+v_vec(2)*gradB_curr(1,2)+v_vec(3)*gradB_curr(1,3))&
                  *dl_norm
        k_1v(2) = (v_vec(1)*gradB_curr(2,1)+v_vec(2)*gradB_curr(2,2)+v_vec(3)*gradB_curr(2,3))&
                  *dl_norm
        k_1v(3) = (v_vec(1)*gradB_curr(3,1)+v_vec(2)*gradB_curr(3,2)+v_vec(3)*gradB_curr(3,3))&
                  *dl_norm

        pos_1in(:) = pos_in(:) + 0.5_num*k_1(:)
        u_1(:) = u_vec(:) + 0.5_num*k_1u(:)
        v_1(:) = v_vec(:) + 0.5_num*k_1v(:)
        call check_position(pos_1in,pos_1out,keep_running)
        call B_gradB_interp(idx_in, pos_1out, B_curr, gradB_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_2(:) = B_curr(:)*dl_norm
        k_2u(1) = (u_1(1)*gradB_curr(1,1)+u_1(2)*gradB_curr(1,2)+u_1(3)*gradB_curr(1,3))*dl_norm
        k_2u(2) = (u_1(1)*gradB_curr(2,1)+u_1(2)*gradB_curr(2,2)+u_1(3)*gradB_curr(2,3))*dl_norm
        k_2u(3) = (u_1(1)*gradB_curr(3,1)+u_1(2)*gradB_curr(3,2)+u_1(3)*gradB_curr(3,3))*dl_norm
        k_2v(1) = (v_1(1)*gradB_curr(1,1)+v_1(2)*gradB_curr(1,2)+v_1(3)*gradB_curr(1,3))*dl_norm
        k_2v(2) = (v_1(1)*gradB_curr(2,1)+v_1(2)*gradB_curr(2,2)+v_1(3)*gradB_curr(2,3))*dl_norm
        k_2v(3) = (v_1(1)*gradB_curr(3,1)+v_1(2)*gradB_curr(3,2)+v_1(3)*gradB_curr(3,3))*dl_norm

        pos_2in(:) = pos_in(:) + 0.5_num*k_2(:)
        u_2(:) = u_vec(:) + 0.5_num*k_2u(:)
        v_2(:) = v_vec(:) + 0.5_num*k_2v(:)
        call check_position(pos_2in,pos_2out,keep_running)
        call B_gradB_interp(idx_in, pos_2out, B_curr, gradB_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_3(:) = B_curr(:)*dl_norm
        k_3u(1) = (u_2(1)*gradB_curr(1,1)+u_2(2)*gradB_curr(1,2)+u_2(3)*gradB_curr(1,3))*dl_norm
        k_3u(2) = (u_2(1)*gradB_curr(2,1)+u_2(2)*gradB_curr(2,2)+u_2(3)*gradB_curr(2,3))*dl_norm
        k_3u(3) = (u_2(1)*gradB_curr(3,1)+u_2(2)*gradB_curr(3,2)+u_2(3)*gradB_curr(3,3))*dl_norm
        k_3v(1) = (v_2(1)*gradB_curr(1,1)+v_2(2)*gradB_curr(1,2)+v_2(3)*gradB_curr(1,3))*dl_norm
        k_3v(2) = (v_2(1)*gradB_curr(2,1)+v_2(2)*gradB_curr(2,2)+v_2(3)*gradB_curr(2,3))*dl_norm
        k_3v(3) = (v_2(1)*gradB_curr(3,1)+v_2(2)*gradB_curr(3,2)+v_2(3)*gradB_curr(3,3))*dl_norm

        pos_3in(:) = pos_in(:) + k_3(:)
        u_3(:) = u_vec(:) + k_3u(:)
        v_3(:) = v_vec(:) + k_3v(:)
        call check_position(pos_3in,pos_3out,keep_running)
        call B_gradB_interp(idx_in, pos_3out, B_curr, gradB_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        k_4(:) = B_curr(:)*dl_norm
        k_4u(1) = (u_3(1)*gradB_curr(1,1)+u_3(2)*gradB_curr(1,2)+u_3(3)*gradB_curr(1,3))*dl_norm
        k_4u(2) = (u_3(1)*gradB_curr(2,1)+u_3(2)*gradB_curr(2,2)+u_3(3)*gradB_curr(2,3))*dl_norm
        k_4u(3) = (u_3(1)*gradB_curr(3,1)+u_3(2)*gradB_curr(3,2)+u_3(3)*gradB_curr(3,3))*dl_norm
        k_4v(1) = (v_3(1)*gradB_curr(1,1)+v_3(2)*gradB_curr(1,2)+v_3(3)*gradB_curr(1,3))*dl_norm
        k_4v(2) = (v_3(1)*gradB_curr(2,1)+v_3(2)*gradB_curr(2,2)+v_3(3)*gradB_curr(2,3))*dl_norm
        k_4v(3) = (v_3(1)*gradB_curr(3,1)+v_3(2)*gradB_curr(3,2)+v_3(3)*gradB_curr(3,3))*dl_norm

        pos_out(:) = pos_in(:) + (k_1(:) + 2.0_num*k_2(:) + 2.0_num*k_3(:) + k_4(:))/6.0_num
        call check_position(pos_out,pos_4out,keep_running)
        IF (keep_running) THEN
          u_vec(:) = u_vec(:) + (k_1u(:) + 2.0_num*k_2u(:) + 2.0_num*k_3u(:) + k_4u(:))/6.0_num
          v_vec(:) = v_vec(:) + (k_1v(:) + 2.0_num*k_2v(:) + 2.0_num*k_3v(:) + k_4v(:))/6.0_num
        END IF

      end subroutine step_cartesianQ_RK4


      subroutine trace_cartesian(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                 B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                 pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace fieldlines, keep only endpoints

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_in(9)
        LOGICAL :: keep_running
        REAL(num) :: pos_out(3), B_out(3), mod_Bout, last_stepsize
        REAL(num) :: pos_next(3), pos_curr(3), u_0(3), v_0(3)

        idx_in(:) = 1

      !Backwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, -dl, &
                               keep_running, check_position, B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                              -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                              B_interp, B_gradB_interp)
        end if
        pos_endpoints(1:3,idx_t) = pos_next(:)

      !Forwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                               dl, keep_running, check_position, B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_in,pos_curr(:), pos_next(:), u_0, v_0, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        pos_endpoints(4:6,idx_t) = pos_next(:)

        if (save_connection) then
          if ((pos_endpoints(3,idx_t) .lt. closed_fl_size)) then                 !one end at photosphere
            if ((pos_endpoints(6,idx_t) .lt. closed_fl_size)) then               !second end at photosphere
              fieldline_connection(idx_t) = 0   !Closed
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 1   !Open
            end if
          else                                                                   !one end away from photosphere
            if ((pos_endpoints(6,idx_t) .lt. closed_fl_size)) then               !second end at photosphere
              fieldline_connection(idx_t) = 1   !Open
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 2   !disconnected
            end if
          end if
        end if

      end subroutine trace_cartesian


      subroutine trace_cartesian_f(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                   B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                   pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace and keep full fieldline locations

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_in(9),step_start,step_total
        LOGICAL :: keep_running
        REAL(num) :: pos_out(3), B_out(3), mod_Bout, last_stepsize, u_0(3), v_0(3)

        pos_fieldline(:,MAX_STEPS+1,idx_t) = pos_start(:)
        idx_in(:) = 1

      !Backwards integration
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          call single_step(idx_in, pos_fieldline(:,MAX_STEPS+2-step_num,idx_t), &
                               pos_fieldline(:,MAX_STEPS+1-step_num,idx_t), u_0, v_0, &
                               -dl, keep_running, check_position, B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        step_start = MAX_STEPS - step_num + 2
        pos_step_start(idx_t) = step_start
        if (.not. keep_running) then
          call check_position(pos_fieldline(:,step_start+1,idx_t),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_fieldline(:,step_start+1,idx_t),B_out,last_stepsize)
          call single_step(idx_in, pos_fieldline(:,step_start+1,idx_t), &
                              pos_fieldline(:,step_start,idx_t), u_0, v_0, &
                              -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                              B_interp, B_gradB_interp)
        end if
        pos_endpoints(1:3,idx_t) = pos_fieldline(:,step_start,idx_t)

      !Forwards integration
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          call single_step(idx_in, pos_fieldline(:,MAX_STEPS+step_num,idx_t), &
                               pos_fieldline(:,MAX_STEPS+1+step_num,idx_t), u_0, v_0, &
                               dl, keep_running, check_position, B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        step_total = MAX_STEPS + step_num - step_start + 1
        pos_step_total(idx_t) = step_total
        if (.not. keep_running) then
          call check_position(pos_fieldline(:,step_total+step_start-2,idx_t),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_fieldline(:,step_total+step_start-2,idx_t),B_out, &
                                  last_stepsize)
          call single_step(idx_in,pos_fieldline(:,step_total+step_start-2,idx_t), &
                               pos_fieldline(:,step_total+step_start-1,idx_t), u_0, v_0, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        pos_endpoints(4:6,idx_t) = pos_fieldline(:,step_total+step_start-1,idx_t)

        if (save_connection) then
          if ((pos_endpoints(3,idx_t) .lt. closed_fl_size)) then                 !one end at photosphere
            if ((pos_endpoints(6,idx_t) .lt. closed_fl_size)) then               !second end at photosphere
              fieldline_connection(idx_t) = 0   !Closed
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 1   !Open
            end if
          else                                                                   !one end away from photosphere
            if ((pos_endpoints(6,idx_t) .lt. closed_fl_size)) then               !second end at photosphere
              fieldline_connection(idx_t) = 1   !Open
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 2   !disconnected
            end if
          end if
        end if

      end subroutine trace_cartesian_f


      subroutine trace_cartesian_Q(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                   B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                   pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace and calculate Q

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_in(9)
        LOGICAL :: keep_running
        REAL(num) :: B_0(3), u_0(3), v_0(3), pos_out(3), B_out(3), mod_Bout, last_stepsize
        REAL(num) :: B_up(3), u_up(3), v_up(3), u_upt(3), v_upt(3)
        REAL(num) :: B_down(3), u_down(3), v_down(3), u_downt(3), v_downt(3)
        REAL(num) :: mod_B0, mod_Bup, mod_Bdown, Q_sign, pos_next(3), pos_curr(3)

        idx_in(:) = 1
        call check_position(pos_start,pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_0)
        mod_B0 = SQRT(B_0(1)**2+B_0(2)**2+B_0(3)**2)
        call create_perpendicular_vectors(B_0,u_0,v_0)

      !Backwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        u_down(:) = u_0(:)
        v_down(:) = v_0(:)
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_down, v_down, &
                               -dl, keep_running, check_position, B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_down, v_down, &
                               -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        pos_endpoints(1:3,idx_t) = pos_next(:)
        call check_position(pos_endpoints(1:3,idx_t),pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_down)
        mod_Bdown = SQRT(B_down(1)**2+B_down(2)**2+B_down(3)**2)
        u_downt(:) = u_down(:) - B_down(:)*vecdot(u_down,B_down)/(mod_Bdown**2)
        v_downt(:) = v_down(:) - B_down(:)*vecdot(v_down,B_down)/(mod_Bdown**2)

      !Forwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        u_up(:) = u_0(:)
        v_up(:) = v_0(:)
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_up, v_up, &
                               dl, keep_running, check_position, B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out, &
                                  last_stepsize)
          call single_step(idx_in,pos_curr(:), pos_next(:), u_down, v_down, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        pos_endpoints(4:6,idx_t) = pos_next(:)
        call check_position(pos_endpoints(4:6,idx_t),pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_up)
        mod_Bup = SQRT(B_up(1)**2+B_up(2)**2+B_up(3)**2)
        u_upt(:) = u_up(:) - B_up(:)*vecdot(u_up,B_up)/(mod_Bup**2)
        v_upt(:) = v_up(:) - B_up(:)*vecdot(v_up,B_up)/(mod_Bup**2)

        if ((pos_endpoints(3,idx_t) .lt. closed_fl_size) .and. &
            (pos_endpoints(6,idx_t) .lt. closed_fl_size)) then
          Q_sign=1.0_num
        else
          Q_sign=-1.0_num
        end if
        !mod-B normalization
        pos_Q(idx_t) = (vecdot(u_upt,u_upt)*vecdot(v_downt,v_downt)+vecdot(v_upt,v_upt)* &
                       vecdot(u_downt,u_downt)-2.0_num*vecdot(u_upt,v_upt) &
                       *vecdot(u_downt,v_downt))*mod_Bup*mod_Bdown/(mod_B0**2)*Q_sign
        !determinant normalization
        !pos_Q(idx_t) = (vecdot(u_upt,u_upt)*vecdot(v_downt,v_downt)+vecdot(v_upt,v_upt)* &
        !               vecdot(u_downt,u_downt)-2.0_num*vecdot(u_upt,v_upt)* &
        !               vecdot(u_downt,v_downt))*Q_sign/SQRT(vecdot(u_downt,u_downt)* &
        !               vecdot(v_downt,v_downt)-(vecdot(u_downt,v_downt))**2)/SQRT( &
        !               vecdot(u_upt,u_upt)*vecdot(v_upt,v_upt)-(vecdot(u_upt,v_upt))**2)

        if (save_connection) then
          if ((pos_endpoints(3,idx_t) .lt. closed_fl_size)) then                 !one end at photosphere
            if ((pos_endpoints(6,idx_t) .lt. closed_fl_size)) then               !second end at photosphere
              fieldline_connection(idx_t) = 0   !Closed
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 1   !Open
            end if
          else                                                                   !one end away from photosphere
            if ((pos_endpoints(6,idx_t) .lt. closed_fl_size)) then               !second end at photosphere
              fieldline_connection(idx_t) = 1   !Open
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 2   !disconnected
            end if
          end if
        end if

      end subroutine trace_cartesian_Q


      subroutine trace_cartesian_Qf(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                    B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                    pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace, calculate Q and keep full fieldline locations

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_in(9),step_start,step_total
        LOGICAL :: keep_running
        REAL(num) :: B_0(3), u_0(3), v_0(3), pos_out(3), B_out(3), mod_Bout, last_stepsize
        REAL(num) :: B_up(3), u_up(3), v_up(3), u_upt(3), v_upt(3)
        REAL(num) :: B_down(3), u_down(3), v_down(3), u_downt(3), v_downt(3)
        REAL(num) :: mod_B0, mod_Bup, mod_Bdown, Q_sign

        pos_fieldline(:,MAX_STEPS+1,idx_t) = pos_start(:)
        idx_in(:) = 1
        call check_position(pos_start,pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_0)
        mod_B0 = SQRT(B_0(1)**2+B_0(2)**2+B_0(3)**2)
        call create_perpendicular_vectors(B_0,u_0,v_0)

      !Backwards integration
        step_num = 1
        keep_running = .true.
        u_down(:) = u_0(:)
        v_down(:) = v_0(:)
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          call single_step(idx_in, pos_fieldline(:,MAX_STEPS+2-step_num,idx_t), &
                               pos_fieldline(:,MAX_STEPS+1-step_num,idx_t), u_down, v_down,  &
                               -dl, keep_running, check_position, B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        step_start = MAX_STEPS - step_num + 2
        pos_step_start(idx_t) = step_start
        if (.not. keep_running) then
          call check_position(pos_fieldline(:,step_start+1,idx_t),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_fieldline(:,step_start+1,idx_t),B_out,last_stepsize)
          call single_step(idx_in, pos_fieldline(:,step_start+1,idx_t), &
                               pos_fieldline(:,step_start,idx_t), u_down, v_down,  &
                               -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        pos_endpoints(1:3,idx_t) = pos_fieldline(:,step_start,idx_t)
        call check_position(pos_endpoints(1:3,idx_t),pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_down)
        mod_Bdown = SQRT(B_down(1)**2+B_down(2)**2+B_down(3)**2)
        u_downt(:) = u_down(:) - B_down(:)*vecdot(u_down,B_down)/(mod_Bdown**2)
        v_downt(:) = v_down(:) - B_down(:)*vecdot(v_down,B_down)/(mod_Bdown**2)

      !Forwards integration
        step_num = 1
        keep_running = .true.
        u_up(:) = u_0(:)
        v_up(:) = v_0(:)
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          call single_step(idx_in, pos_fieldline(:,MAX_STEPS+step_num,idx_t), &
                               pos_fieldline(:,MAX_STEPS+1+step_num,idx_t), u_up, v_up,   &
                               dl, keep_running, check_position,B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        step_total = MAX_STEPS + step_num - step_start + 1
        pos_step_total(idx_t) = step_total
        if (.not. keep_running) then
          call check_position(pos_fieldline(:,step_total+step_start-2,idx_t),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_fieldline(:,step_total+step_start-2,idx_t),B_out, &
                                  last_stepsize)
          call single_step(idx_in,pos_fieldline(:,step_total+step_start-2,idx_t), &
                               pos_fieldline(:,step_total+step_start-1,idx_t), u_down, v_down,  &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        pos_endpoints(4:6,idx_t) = pos_fieldline(:,step_total+step_start-1,idx_t)
        call check_position(pos_endpoints(4:6,idx_t),pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_up)
        mod_Bup = SQRT(B_up(1)**2+B_up(2)**2+B_up(3)**2)
        u_upt(:) = u_up(:) - B_up(:)*vecdot(u_up,B_up)/(mod_Bup**2)
        v_upt(:) = v_up(:) - B_up(:)*vecdot(v_up,B_up)/(mod_Bup**2)

        if ((pos_endpoints(3,idx_t) .lt. closed_fl_size) .and. &
            (pos_endpoints(6,idx_t).lt. closed_fl_size)) then
          Q_sign=1.0_num
        else
          Q_sign=-1.0_num
        end if
        !mod-B normalization
        pos_Q(idx_t) = (vecdot(u_upt,u_upt)*vecdot(v_downt,v_downt)+vecdot(v_upt,v_upt)* &
                       vecdot(u_downt,u_downt)-2.0_num*vecdot(u_upt,v_upt) &
                       *vecdot(u_downt,v_downt))*mod_Bup*mod_Bdown/(mod_B0**2)*Q_sign
        !determinant normalization
        !pos_Q(idx_t) = (vecdot(u_upt,u_upt)*vecdot(v_downt,v_downt)+vecdot(v_upt,v_upt)* &
        !               vecdot(u_downt,u_downt)-2.0_num*vecdot(u_upt,v_upt)* &
        !               vecdot(u_downt,v_downt))*Q_sign/SQRT(vecdot(u_downt,u_downt)* &
        !               vecdot(v_downt,v_downt)-(vecdot(u_downt,v_downt))**2)/SQRT( &
        !               vecdot(u_upt,u_upt)*vecdot(v_upt,v_upt)-(vecdot(u_upt,v_upt))**2)

        if (save_connection) then
          if ((pos_endpoints(3,idx_t) .lt. closed_fl_size)) then                 !one end at photosphere
            if ((pos_endpoints(6,idx_t) .lt. closed_fl_size)) then               !second end at photosphere
              fieldline_connection(idx_t) = 0   !Closed
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 1   !Open
            end if
          else                                                                   !one end away from photosphere
            if ((pos_endpoints(6,idx_t) .lt. closed_fl_size)) then               !second end at photosphere
              fieldline_connection(idx_t) = 1   !Open
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 2   !disconnected
            end if
          end if
        end if

      end subroutine trace_cartesian_Qf


      subroutine step_spherical(idx_in, pos_in, pos_out, u_vec, v_vec, dl, &
                                keep_running, check_position, B_interp, B_gradB_interp)

        INTEGER :: idx_in(9)
        REAL(num) :: pos_in(3),pos_out(3),u_vec(3),v_vec(3),dl
        LOGICAL :: keep_running
        procedure(check_position_iface) :: check_position
        procedure(B_interp_iface) :: B_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp

        REAL(num) :: mod_B, pos_0out(3), k_1(3), B_curr(3), dl_norm, r_sin_th

        call check_position(pos_in,pos_0out,keep_running)
        call B_interp(idx_in, pos_0out, B_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        r_sin_th = pos_0out(1)*SIN(pos_0out(2))
        k_1(1) = B_curr(1)*dl_norm
        k_1(2) = B_curr(2)*dl_norm/pos_0out(1)
        k_1(3) = B_curr(3)*dl_norm/r_sin_th

        pos_out(:) = pos_in(:) + k_1(:)
        call check_position(pos_out,pos_0out,keep_running)

      end subroutine step_spherical


      subroutine step_sphericalQnoc(idx_in, pos_in, pos_out, u_vec, v_vec, dl, &
                                 keep_running, check_position, B_interp, B_gradB_interp)
      !This subroutine is now deprecated; in principle, it should approach QSLsquasher
      !Note - gradB_curr is just partial derivatives, not actual gradient in sphericals

        INTEGER :: idx_in(9)
        REAL(num) :: pos_in(3),pos_out(3),u_vec(3),v_vec(3),dl
        LOGICAL :: keep_running
        procedure(check_position_iface) :: check_position
        procedure(B_interp_iface) :: B_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp

        REAL(num) :: mod_B, pos_0out(3), k_1(3), B_curr(3), dl_norm
        REAL(num) :: k_1u(3), k_1v(3), gradB_curr(3,3), r_sin_th, cot_th

        call check_position(pos_in,pos_0out,keep_running)
        call B_gradB_interp(idx_in, pos_0out, B_curr, gradB_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        r_sin_th = pos_0out(1)*SIN(pos_0out(2))
        k_1(1) = B_curr(1)*dl_norm
        k_1(2) = B_curr(2)*dl_norm/pos_0out(1)
        k_1(3) = B_curr(3)*dl_norm/r_sin_th
      !Directional derivative as defined: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
        !cot_th = COS(pos_0out(2))/SIN(pos_0out(2))
        !k_1u(1) = (u_vec(1)*gradB_curr(1,1)+u_vec(2)*gradB_curr(1,2)/pos_0out(1) &
        !          +u_vec(3)*gradB_curr(1,3)/r_sin_th-(u_vec(2)*B_curr(2)+u_vec(3)*B_curr(3)) &
        !          /pos_0out(1))*dl_norm
        !k_1u(2) = (u_vec(1)*gradB_curr(2,1)+u_vec(2)*gradB_curr(2,2)/pos_0out(1) &
        !          +u_vec(3)*gradB_curr(2,3)/r_sin_th+(u_vec(2)*B_curr(1)-u_vec(3)*B_curr(3) &
        !          *cot_th)/pos_0out(1))*dl_norm/pos_0out(1)
        !k_1u(3) = (u_vec(1)*gradB_curr(3,1)+u_vec(2)*gradB_curr(3,2)/pos_0out(1) &
        !          +u_vec(3)*gradB_curr(3,3)/r_sin_th+(u_vec(3)*B_curr(1)-u_vec(3)*B_curr(2) &
        !          *cot_th)/pos_0out(1))*dl_norm/r_sin_th
        !k_1v(1) = (v_vec(1)*gradB_curr(1,1)+v_vec(2)*gradB_curr(1,2)/pos_0out(1) &
        !          +v_vec(3)*gradB_curr(1,3)/r_sin_th-(v_vec(2)*B_curr(2)+v_vec(3)*B_curr(3)) &
        !          /pos_0out(1))*dl_norm
        !k_1v(2) = (v_vec(1)*gradB_curr(2,1)+v_vec(2)*gradB_curr(2,2)/pos_0out(1) &
        !          +v_vec(3)*gradB_curr(2,3)/r_sin_th+(v_vec(2)*B_curr(1)-v_vec(3)*B_curr(3) &
        !          *cot_th)/pos_0out(1))*dl_norm/pos_0out(1)
        !k_1v(3) = (v_vec(1)*gradB_curr(3,1)+v_vec(2)*gradB_curr(3,2)/pos_0out(1) &
        !          +v_vec(3)*gradB_curr(3,3)/r_sin_th+(v_vec(3)*B_curr(1)-v_vec(3)*B_curr(2) &
        !          *cot_th)/pos_0out(1))*dl_norm/r_sin_th
      !Surprisingly, this works better than the directional derivative
        k_1u(1) = (u_vec(1)*gradB_curr(1,1)+u_vec(2)*gradB_curr(1,2)+u_vec(3)*gradB_curr(1,3))&
                  *dl_norm
        k_1u(2) = (u_vec(1)*gradB_curr(2,1)+u_vec(2)*gradB_curr(2,2)+u_vec(3)*gradB_curr(2,3))&
                  *dl_norm/pos_0out(1)
        k_1u(3) = (u_vec(1)*gradB_curr(3,1)+u_vec(2)*gradB_curr(3,2)+u_vec(3)*gradB_curr(3,3))&
                  *dl_norm/r_sin_th
        k_1v(1) = (v_vec(1)*gradB_curr(1,1)+v_vec(2)*gradB_curr(1,2)+v_vec(3)*gradB_curr(1,3))&
                  *dl_norm
        k_1v(2) = (v_vec(1)*gradB_curr(2,1)+v_vec(2)*gradB_curr(2,2)+v_vec(3)*gradB_curr(2,3))&
                  *dl_norm/pos_0out(1)
        k_1v(3) = (v_vec(1)*gradB_curr(3,1)+v_vec(2)*gradB_curr(3,2)+v_vec(3)*gradB_curr(3,3))&
                  *dl_norm/r_sin_th
      !Another approach to taking the gradient
        !k_1u(1) = (u_vec(1)*gradB_curr(1,1)+u_vec(2)*gradB_curr(1,2)/pos_0out(1) &
        !          +u_vec(3)*gradB_curr(1,3)/r_sin_th)*dl_norm
        !k_1u(2) = (u_vec(1)*gradB_curr(2,1)+u_vec(2)*gradB_curr(2,2)/pos_0out(1) &
        !          +u_vec(3)*gradB_curr(2,3)/r_sin_th)*dl_norm
        !k_1u(3) = (u_vec(1)*gradB_curr(3,1)+u_vec(2)*gradB_curr(3,2)/pos_0out(1) &
        !          +u_vec(3)*gradB_curr(3,3)/r_sin_th)*dl_norm
        !k_1v(1) = (v_vec(1)*gradB_curr(1,1)+v_vec(2)*gradB_curr(1,2)/pos_0out(1) &
        !          +v_vec(3)*gradB_curr(1,3)/r_sin_th)*dl_norm
        !k_1v(2) = (v_vec(1)*gradB_curr(2,1)+v_vec(2)*gradB_curr(2,2)/pos_0out(1) &
        !          +v_vec(3)*gradB_curr(2,3)/r_sin_th)*dl_norm
        !k_1v(3) = (v_vec(1)*gradB_curr(3,1)+v_vec(2)*gradB_curr(3,2)/pos_0out(1) &
        !          +v_vec(3)*gradB_curr(3,3)/r_sin_th)*dl_norm

        pos_out(:) = pos_in(:) + k_1(:)
        call check_position(pos_out,pos_0out,keep_running)
        IF (keep_running) THEN
          u_vec(:) = u_vec(:) + k_1u(:)
          v_vec(:) = v_vec(:) + k_1v(:)
        END IF

      end subroutine step_sphericalQnoc


      subroutine step_sphericalQ(idx_in, pos_in, pos_out, u_vec, v_vec, dl, &
                                 keep_running, check_position, B_interp, B_gradB_interp)
      !Note - gradB_curr is just partial derivatives, not actual gradient in sphericals
      !include curvature 

        INTEGER :: idx_in(9)
        REAL(num) :: pos_in(3),pos_out(3),u_vec(3),v_vec(3),dl
        LOGICAL :: keep_running
        procedure(check_position_iface) :: check_position
        procedure(B_interp_iface) :: B_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp

        REAL(num) :: mod_B, pos_0out(3), k_1(3), B_curr(3), dl_norm
        REAL(num) :: k_1u(3), k_1v(3), gradB_curr(3,3), r_sin_th, cot_th

        call check_position(pos_in,pos_0out,keep_running)
        call B_gradB_interp(idx_in, pos_0out, B_curr, gradB_curr)
        mod_B = SQRT(B_curr(1)**2+B_curr(2)**2+B_curr(3)**2)
        dl_norm = dl/mod_B !Step size scaled by B
        r_sin_th = pos_0out(1)*SIN(pos_0out(2))
        k_1(1) = B_curr(1)*dl_norm
        k_1(2) = B_curr(2)*dl_norm/pos_0out(1)
        k_1(3) = B_curr(3)*dl_norm/r_sin_th
      !Formula with curvature
        cot_th = COS(pos_0out(2))/SIN(pos_0out(2))
        k_1u(1) = (u_vec(1)*gradB_curr(1,1)+u_vec(2)*gradB_curr(1,2)/pos_0out(1) &
                  +u_vec(3)*gradB_curr(1,3)/r_sin_th)*dl_norm
        k_1u(2) = (u_vec(1)*gradB_curr(2,1)+u_vec(2)*gradB_curr(2,2)/pos_0out(1) &
                  +u_vec(3)*gradB_curr(2,3)/r_sin_th+(u_vec(2)*B_curr(1)-u_vec(1)*B_curr(2) &
                  )/pos_0out(1))*dl_norm
        k_1u(3) = (u_vec(1)*gradB_curr(3,1)+u_vec(2)*gradB_curr(3,2)/pos_0out(1) &
                  +u_vec(3)*gradB_curr(3,3)/r_sin_th+(u_vec(3)*B_curr(1)-u_vec(1)*B_curr(3) &
                  +(u_vec(3)*B_curr(2)-u_vec(2)*B_curr(3))*cot_th)/pos_0out(1))*dl_norm
        k_1v(1) = (v_vec(1)*gradB_curr(1,1)+v_vec(2)*gradB_curr(1,2)/pos_0out(1) &
                  +v_vec(3)*gradB_curr(1,3)/r_sin_th)*dl_norm
        k_1v(2) = (v_vec(1)*gradB_curr(2,1)+v_vec(2)*gradB_curr(2,2)/pos_0out(1) &
                  +v_vec(3)*gradB_curr(2,3)/r_sin_th+(v_vec(2)*B_curr(1)-v_vec(1)*B_curr(2) &
                  )/pos_0out(1))*dl_norm
        k_1v(3) = (v_vec(1)*gradB_curr(3,1)+v_vec(2)*gradB_curr(3,2)/pos_0out(1) &
                  +v_vec(3)*gradB_curr(3,3)/r_sin_th+(v_vec(3)*B_curr(1)-v_vec(1)*B_curr(3) &
                  +(v_vec(3)*B_curr(2)-v_vec(2)*B_curr(3))*cot_th)/pos_0out(1))*dl_norm

        pos_out(:) = pos_in(:) + k_1(:)
        call check_position(pos_out,pos_0out,keep_running)
        IF (keep_running) THEN
          u_vec(:) = u_vec(:) + k_1u(:)
          v_vec(:) = v_vec(:) + k_1v(:)
        END IF

      end subroutine step_sphericalQ


      subroutine trace_spherical(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                 B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                 pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace fieldlines, keep only endpoints

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_in(9)
        LOGICAL :: keep_running
        REAL(num) :: pos_out(3), B_out(3), mod_Bout, last_stepsize
        REAL(num) :: pos_next(3), pos_curr(3), u_0(3), v_0(3)

        idx_in(:) = 1

      !Backwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, -dl, &
                               keep_running, check_position,B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                              -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                              B_interp, B_gradB_interp)
        end if
        pos_endpoints(1:3,idx_t) = pos_next(:)

      !Forwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                               dl, keep_running, check_position,B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        pos_endpoints(4:6,idx_t) = pos_next(:)

        if (save_connection) then
          if ((pos_endpoints(1,idx_t) .lt. closed_fl_size)) then                 !one end at photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              fieldline_connection(idx_t) = 0   !Closed
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 1   !Open
            end if
          else                                                                   !one end away from photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              fieldline_connection(idx_t) = 1   !Open
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 2   !disconnected
            end if
          end if
        end if

      end subroutine trace_spherical


      subroutine trace_spherical_f(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                   B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                   pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace and keep full fieldline locations

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_in(9), step_start, step_total
        LOGICAL :: keep_running
        REAL(num) :: pos_out(3), B_out(3), mod_Bout, last_stepsize, u_0(3), v_0(3)
        REAL(num) :: pos_next(3), pos_curr(3)

        pos_fieldline(:,MAX_STEPS+1,idx_t) = pos_start(:)
        idx_in(:) = 1

      !Backwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, -dl, &
                               keep_running, check_position,B_interp, B_gradB_interp)
          pos_fieldline(:,MAX_STEPS+1-step_num,idx_t) = pos_next(:)
          step_num = step_num + 1
        end do
        step_start = MAX_STEPS - step_num + 2
        pos_step_start(idx_t) = step_start
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                              -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                              B_interp, B_gradB_interp)
          pos_fieldline(:,step_start,idx_t) = pos_next(:)
        end if
        pos_endpoints(1:3,idx_t) = pos_next(:)

      !Forwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                               dl, keep_running, check_position,B_interp, B_gradB_interp)
          pos_fieldline(:,MAX_STEPS+1+step_num,idx_t) = pos_next(:)
          step_num = step_num + 1
        end do
        step_total = MAX_STEPS + step_num - step_start + 1
        pos_step_total(idx_t) = step_total
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
          pos_fieldline(:,step_total+step_start-1,idx_t) = pos_next(:)
        end if
        pos_endpoints(4:6,idx_t) = pos_next(:)

        if (save_connection) then
          if ((pos_endpoints(1,idx_t) .lt. closed_fl_size)) then                 !one end at photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              fieldline_connection(idx_t) = 0   !Closed
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 1   !Open
            end if
          else                                                                   !one end away from photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              fieldline_connection(idx_t) = 1   !Open
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 2   !disconnected
            end if
          end if
        end if

      end subroutine trace_spherical_f


      subroutine trace_spherical_Q(check_position,intercept_boundary,B_interp,Bfull_interp,&
                                   B_gradB_interp, single_step,pos_start,idx_t,pos_endpoints, &
                                   pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace and calculate Q

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_in(9)
        LOGICAL :: keep_running
        REAL(num) :: B_0(3), u_0(3), v_0(3), pos_out(3), B_out(3), mod_Bout, last_stepsize
        REAL(num) :: B_up(3), u_up(3), v_up(3), u_upt(3), v_upt(3)
        REAL(num) :: B_down(3), u_down(3), v_down(3), u_downt(3), v_downt(3)
        REAL(num) :: mod_B0, mod_Bup, mod_Bdown, Q_sign, pos_next(3), pos_curr(3)

        idx_in(:) = 1
        call check_position(pos_start,pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_0)
        mod_B0 = SQRT(B_0(1)**2+B_0(2)**2+B_0(3)**2)
        call create_perpendicular_vectors(B_0,u_0,v_0)

      !Backwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        u_down(:) = u_0(:)
        v_down(:) = v_0(:)
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_down, v_down, &
                               -dl, keep_running, check_position,B_interp,B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_down, v_down, &
                               -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp,B_gradB_interp)
        end if
        pos_endpoints(1:3,idx_t) = pos_next(:)
        call check_position(pos_endpoints(1:3,idx_t),pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_down)
        mod_Bdown = SQRT(B_down(1)**2+B_down(2)**2+B_down(3)**2)
        u_downt(:) = u_down(:) - B_down(:)*vecdot(u_down,B_down)/(mod_Bdown**2)
        v_downt(:) = v_down(:) - B_down(:)*vecdot(v_down,B_down)/(mod_Bdown**2)

      !Forwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        u_up(:) = u_0(:)
        v_up(:) = v_0(:)
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_up, v_up, &
                               dl, keep_running, check_position,B_interp,B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out, &
                                  last_stepsize)
          call single_step(idx_in,pos_curr(:), pos_next(:), u_up, v_up, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp,B_gradB_interp)
        end if
        pos_endpoints(4:6,idx_t) = pos_next(:)
        call check_position(pos_endpoints(4:6,idx_t),pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_up)
        mod_Bup = SQRT(B_up(1)**2+B_up(2)**2+B_up(3)**2)
        u_upt(:) = u_up(:) - B_up(:)*vecdot(u_up,B_up)/(mod_Bup**2)
        v_upt(:) = v_up(:) - B_up(:)*vecdot(v_up,B_up)/(mod_Bup**2)

        if ((pos_endpoints(1,idx_t) .lt. closed_fl_size) .and. &
            (pos_endpoints(4,idx_t).lt. closed_fl_size)) then
          Q_sign=1.0_num
        else
          Q_sign=-1.0_num
        end if
        !mod-B normalization
        pos_Q(idx_t) = (vecdot(u_upt,u_upt)*vecdot(v_downt,v_downt)+vecdot(v_upt,v_upt)* &
                        vecdot(u_downt,u_downt)-2.0_num*vecdot(u_upt,v_upt) &
                        *vecdot(u_downt,v_downt))*Q_sign*mod_Bup*mod_Bdown/(mod_B0**2)
        !determinant normalization
        !pos_Q(idx_t) = (vecdot(u_upt,u_upt)*vecdot(v_downt,v_downt)+vecdot(v_upt,v_upt)* &
        !               vecdot(u_downt,u_downt)-2.0_num*vecdot(u_upt,v_upt)* &
        !               vecdot(u_downt,v_downt))*Q_sign/SQRT(vecdot(u_downt,u_downt)* &
        !               vecdot(v_downt,v_downt)-(vecdot(u_downt,v_downt))**2)/SQRT( &
        !               vecdot(u_upt,u_upt)*vecdot(v_upt,v_upt)-(vecdot(u_upt,v_upt))**2)

        if (save_connection) then
          if ((pos_endpoints(1,idx_t) .lt. closed_fl_size)) then                 !one end at photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              fieldline_connection(idx_t) = 0   !Closed
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 1   !Open
            end if
          else                                                                   !one end away from photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              fieldline_connection(idx_t) = 1   !Open
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 2   !disconnected
            end if
          end if
        end if

      end subroutine trace_spherical_Q


      subroutine trace_spherical_Qf(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                    B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                    pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace, calculate Q and keep full fieldline locations

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_in(9), step_start, step_total
        LOGICAL :: keep_running
        REAL(num) :: B_0(3), u_0(3), v_0(3), pos_out(3), B_out(3), mod_Bout, last_stepsize
        REAL(num) :: B_up(3), u_up(3), v_up(3), u_upt(3), v_upt(3)
        REAL(num) :: B_down(3), u_down(3), v_down(3), u_downt(3), v_downt(3)
        REAL(num) :: mod_B0, mod_Bup, mod_Bdown, Q_sign, pos_next(3), pos_curr(3)

        pos_fieldline(:,MAX_STEPS+1,idx_t) = pos_start(:)
        idx_in(:) = 1
        call check_position(pos_start,pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_0)
        mod_B0 = SQRT(B_0(1)**2+B_0(2)**2+B_0(3)**2)
        call create_perpendicular_vectors(B_0,u_0,v_0)

      !Backwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        u_down(:) = u_0(:)
        v_down(:) = v_0(:)
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_down, v_down, &
                               -dl, keep_running, check_position,B_interp,B_gradB_interp)
          pos_fieldline(:,MAX_STEPS+1-step_num,idx_t) = pos_next(:)
          step_num = step_num + 1
        end do
        step_start = MAX_STEPS - step_num + 2
        pos_step_start(idx_t) = step_start
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_down, v_down, &
                               -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp,B_gradB_interp)
          pos_fieldline(:,step_start,idx_t) = pos_next(:)
        end if
        pos_endpoints(1:3,idx_t) = pos_next(:)
        call check_position(pos_endpoints(1:3,idx_t),pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_down)
        mod_Bdown = SQRT(B_down(1)**2+B_down(2)**2+B_down(3)**2)
        u_downt(:) = u_down(:) - B_down(:)*vecdot(u_down,B_down)/(mod_Bdown**2)
        v_downt(:) = v_down(:) - B_down(:)*vecdot(v_down,B_down)/(mod_Bdown**2)

      !Forwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        u_up(:) = u_0(:)
        v_up(:) = v_0(:)
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_in, pos_curr(:), pos_next(:), u_up, v_up, &
                               dl, keep_running, check_position,B_interp,B_gradB_interp)
          pos_fieldline(:,MAX_STEPS+1+step_num,idx_t) = pos_next(:)
          step_num = step_num + 1
        end do
        step_total = MAX_STEPS + step_num - step_start + 1
        pos_step_total(idx_t) = step_total
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_in, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out, &
                                  last_stepsize)
          call single_step(idx_in,pos_curr(:), pos_next(:), u_up, v_up, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp,B_gradB_interp)
          pos_fieldline(:,step_total+step_start-1,idx_t) = pos_next(:)
        end if
        pos_endpoints(4:6,idx_t) = pos_next(:)
        call check_position(pos_endpoints(4:6,idx_t),pos_out,keep_running)
        call Bfull_interp(idx_in, pos_out, B_up)
        mod_Bup = SQRT(B_up(1)**2+B_up(2)**2+B_up(3)**2)
        u_upt(:) = u_up(:) - B_up(:)*vecdot(u_up,B_up)/(mod_Bup**2)
        v_upt(:) = v_up(:) - B_up(:)*vecdot(v_up,B_up)/(mod_Bup**2)

        if ((pos_endpoints(1,idx_t) .lt. closed_fl_size) .and. &
            (pos_endpoints(4,idx_t).lt. closed_fl_size)) then
          Q_sign=1.0_num
        else
          Q_sign=-1.0_num
        end if
        !mod-B normalization
        pos_Q(idx_t) = (vecdot(u_upt,u_upt)*vecdot(v_downt,v_downt)+vecdot(v_upt,v_upt)* &
                        vecdot(u_downt,u_downt)-2.0_num*vecdot(u_upt,v_upt) &
                        *vecdot(u_downt,v_downt))*Q_sign*mod_Bup*mod_Bdown/(mod_B0**2)
        !determinant normalization
        !pos_Q(idx_t) = (vecdot(u_upt,u_upt)*vecdot(v_downt,v_downt)+vecdot(v_upt,v_upt)* &
        !               vecdot(u_downt,u_downt)-2.0_num*vecdot(u_upt,v_upt)* &
        !               vecdot(u_downt,v_downt))*Q_sign/SQRT(vecdot(u_downt,u_downt)* &
        !               vecdot(v_downt,v_downt)-(vecdot(u_downt,v_downt))**2)/SQRT( &
        !               vecdot(u_upt,u_upt)*vecdot(v_upt,v_upt)-(vecdot(u_upt,v_upt))**2)

        if (save_connection) then
          if ((pos_endpoints(1,idx_t) .lt. closed_fl_size)) then                 !one end at photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              fieldline_connection(idx_t) = 0   !Closed
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 1   !Open
            end if
          else                                                                   !one end away from photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              fieldline_connection(idx_t) = 1   !Open
            else                                                                 !second end away from photosphere
              fieldline_connection(idx_t) = 2   !disconnected
            end if
          end if
        end if

      end subroutine trace_spherical_Qf


      subroutine run_trace

        INTEGER :: idx_t
        REAL(num) :: fieldline_start(3)

        procedure (get_position0_iface), pointer :: gpos_ptr => null ()
        procedure (check_position_iface), pointer :: cpos_ptr => null ()
        procedure (intercept_boundary_iface), pointer :: ibdry_ptr => null ()
        !B_interp can use B_normalized (before interpolation), or not
        procedure (B_interp_iface), pointer :: binterp_ptr => null ()
        !Bfull_interp must use B unnormalized
        procedure (Bfull_interp_iface), pointer :: bfullinterp_ptr => null ()
        procedure (B_gradB_interp_iface), pointer :: bgradinterp_ptr => null ()
        procedure (single_step_iface), pointer :: step_ptr => null ()
        procedure (trace_fl_iface), pointer :: tr_ptr => null ()

        if (user_defined) then
          save_endpoints = .false.
          save_Q = .false.
          save_connection = .false.
          save_fieldlines = .false.
        end if

        if ((.not. save_endpoints) .and. (.not. save_Q) .and. (.not. save_connection) .and. &
            (.not. save_fieldlines) .and. (.not. user_defined)) then
          save_endpoints = .true.
        end if

        ALLOCATE(fieldline_endpoints(6,numin_tot))
        if (save_fieldlines) then
          ALLOCATE(fieldline_allpos(3,2*MAX_STEPS+1,numin_tot))
          ALLOCATE(fieldline_pts(numin_tot))
          ALLOCATE(fieldline_ptn(numin_tot))
        end if
        if (save_Q) then
          ALLOCATE(fieldline_Q(numin_tot))
        end if
        if (save_connection) then
          ALLOCATE(fieldline_connection(numin_tot))
        end if
        if (user_defined) then
          call prepare_user_defined
        end if

        if (input_type .eq. 0) then
          gpos_ptr => get_start_pos_0
        else if ((input_type .eq. 1) .or. (input_type .eq. 2)) then
          gpos_ptr => get_start_pos_12
        end if

        if (grid_separate) then
          if (normalized_B) then
            print *, 'Warning: normalized B (-nb) is incompatible with'
            print *, 'grid separate (-gs); field will not be normalized'
          end if
          if (grid_regular) then
            binterp_ptr => B_interp_regular_separate
            bfullinterp_ptr => B_interp_regular_separate
            bgradinterp_ptr => B_gradB_interp_regular_separate
          else
            binterp_ptr => B_interp_irregular
            bfullinterp_ptr => B_interp_irregular
            bgradinterp_ptr => B_gradB_interp_irregular
          end if
        else
          if ((grid_regular) .and. (normalized_B)) then
            binterp_ptr => B_interp_normalized_regular
            bfullinterp_ptr => B_interp_regular
            bgradinterp_ptr => B_gradB_interp_normalized_regular
          else if ((grid_regular) .and. (.not. normalized_B)) then
            binterp_ptr => B_interp_regular
            bfullinterp_ptr => B_interp_regular
            bgradinterp_ptr => B_gradB_interp_regular
          else if ((.not. grid_regular) .and. (normalized_B)) then
            binterp_ptr => B_interp_normalized_irregular
            bfullinterp_ptr => B_interp_irregular
            bgradinterp_ptr => B_gradB_interp_normalized_irregular
          else
            binterp_ptr => B_interp_irregular
            bfullinterp_ptr => B_interp_irregular
            bgradinterp_ptr => B_gradB_interp_irregular
          end if
        end if

        if (geometry .eq. 0) then !Cartesian

          if (user_defined) then
            step_ptr => step_cartesian
            tr_ptr => trace_cartesian_user
          else if (save_Q .and. save_fieldlines) then
            step_ptr => step_cartesianQ
            tr_ptr => trace_cartesian_Qf
          else if (save_Q .and. (.not. save_fieldlines)) then
            step_ptr => step_cartesianQ
            tr_ptr => trace_cartesian_Q
          else if ((.not. save_Q) .and. save_fieldlines) then
            step_ptr => step_cartesian
            tr_ptr => trace_cartesian_f
          else if ((.not. save_Q) .and. (.not. save_fieldlines)) then
            step_ptr => step_cartesian
            tr_ptr => trace_cartesian
          end if

          if (periodic_X .and. periodic_Y .and. periodic_Z) then
            cpos_ptr => check_position_c111
            ibdry_ptr => intercept_boundary_c111
          else if (periodic_X .and. periodic_Y .and. (.not. periodic_Z)) then
            cpos_ptr => check_position_c110
            ibdry_ptr => intercept_boundary_c110
          else if (periodic_X .and. (.not. periodic_Y) .and. periodic_Z) then
            cpos_ptr => check_position_c101
            ibdry_ptr => intercept_boundary_c101
          else if (periodic_X .and. (.not. periodic_Y) .and. (.not. periodic_Z)) then
            cpos_ptr => check_position_c100
            ibdry_ptr => intercept_boundary_c100
          else if ((.not. periodic_X) .and. periodic_Y .and. periodic_Z) then
            cpos_ptr => check_position_c011
            ibdry_ptr => intercept_boundary_c011
          else if ((.not. periodic_X) .and. periodic_Y .and. (.not. periodic_Z)) then
            cpos_ptr => check_position_c010
            ibdry_ptr => intercept_boundary_c010
          else if ((.not. periodic_X) .and. (.not. periodic_Y) .and. periodic_Z) then
            cpos_ptr => check_position_c001
            ibdry_ptr => intercept_boundary_c001
          else if ((.not. periodic_X) .and. (.not. periodic_Y) .and. (.not. periodic_Z)) then
            cpos_ptr => check_position_c000
            ibdry_ptr => intercept_boundary_c000
          end if

          closed_fl_size=closed_fl_constant*(grid3max-grid3min)+grid3min

        else if (geometry .eq. 1) then !Spherical

          if (user_defined) then
            step_ptr => step_spherical
            tr_ptr => trace_spherical_user
          else if (save_Q .and. save_fieldlines) then
            step_ptr => step_sphericalQ
            tr_ptr => trace_spherical_Qf
          else if (save_Q .and. (.not. save_fieldlines)) then
            step_ptr => step_sphericalQ
            tr_ptr => trace_spherical_Q
          else if ((.not. save_Q) .and. save_fieldlines) then
            step_ptr => step_spherical
            tr_ptr => trace_spherical_f
          else if ((.not. save_Q) .and. (.not. save_fieldlines)) then
            step_ptr => step_spherical
            tr_ptr => trace_spherical
          end if

          !If theta boundary is within theta_periodicity_threshold, it's open - field lines can pass over pole if necessary
          if (periodic_PHI) then
            if ((grid2min .lt. theta_periodicity_threshold) .and. &
                (grid2max .gt. PI-theta_periodicity_threshold)) then
              cpos_ptr => check_position_s031
              ibdry_ptr => intercept_boundary_s031
            else if (grid2min .lt. theta_periodicity_threshold) then
              cpos_ptr => check_position_s011
              ibdry_ptr => intercept_boundary_s011
            else if (grid2max .gt. PI-theta_periodicity_threshold) then
              cpos_ptr => check_position_s021
              ibdry_ptr => intercept_boundary_s021
            else
              cpos_ptr => check_position_s001
              ibdry_ptr => intercept_boundary_s001
            end if
          else if ((.not. periodic_PHI)) then
            if ((grid2min .lt. theta_periodicity_threshold) .and. &
                (grid2max .gt. PI-theta_periodicity_threshold)) then
              cpos_ptr => check_position_s030
              ibdry_ptr => intercept_boundary_s030
            else if (grid2min .lt. theta_periodicity_threshold) then
              cpos_ptr => check_position_s010
              ibdry_ptr => intercept_boundary_s010
            else if (grid2max .gt. PI-theta_periodicity_threshold) then
              cpos_ptr => check_position_s020
              ibdry_ptr => intercept_boundary_s020
            else
              cpos_ptr => check_position_s000
              ibdry_ptr => intercept_boundary_s000
            end if
          end if

          closed_fl_size=closed_fl_constant*(grid1max-grid1min)+grid1min

        end if

      !Main tracing loop
      !$OMP parallel num_threads(num_proc) private(fieldline_start)
      !$OMP do
        do idx_t = 1, numin_tot
          call gpos_ptr(fieldline_start,idx_t)
          call tr_ptr(cpos_ptr,ibdry_ptr,binterp_ptr,bfullinterp_ptr,bgradinterp_ptr,step_ptr, &
                      fieldline_start,idx_t,fieldline_endpoints, fieldline_Q, fieldline_allpos, &
                      fieldline_pts, fieldline_ptn, step_size)
        end do
      !$OMP end do
      !$OMP end parallel

      end subroutine run_trace


      subroutine cleanup

        IF (ALLOCATED(coord1_in)) DEALLOCATE(coord1_in)
        IF (ALLOCATED(coord2_in)) DEALLOCATE(coord2_in)
        IF (ALLOCATED(coord3_in)) DEALLOCATE(coord3_in)
        IF (ALLOCATED(grid1)) DEALLOCATE(grid1)
        IF (ALLOCATED(grid2)) DEALLOCATE(grid2)
        IF (ALLOCATED(grid3)) DEALLOCATE(grid3)
        IF (ALLOCATED(B_grid)) DEALLOCATE(B_grid)
        IF (ALLOCATED(grid1_1)) DEALLOCATE(grid1_1)
        IF (ALLOCATED(grid1_2)) DEALLOCATE(grid1_2)
        IF (ALLOCATED(grid1_3)) DEALLOCATE(grid1_3)
        IF (ALLOCATED(grid2_1)) DEALLOCATE(grid2_1)
        IF (ALLOCATED(grid2_2)) DEALLOCATE(grid2_2)
        IF (ALLOCATED(grid2_3)) DEALLOCATE(grid2_3)
        IF (ALLOCATED(grid3_1)) DEALLOCATE(grid3_1)
        IF (ALLOCATED(grid3_2)) DEALLOCATE(grid3_2)
        IF (ALLOCATED(grid3_3)) DEALLOCATE(grid3_3)
        IF (ALLOCATED(B_grid1)) DEALLOCATE(B_grid1)
        IF (ALLOCATED(B_grid2)) DEALLOCATE(B_grid2)
        IF (ALLOCATED(B_grid3)) DEALLOCATE(B_grid3)
        IF (ALLOCATED(grid1_ir)) DEALLOCATE(grid1_ir)
        IF (ALLOCATED(grid2_ir)) DEALLOCATE(grid2_ir)
        IF (ALLOCATED(grid3_ir)) DEALLOCATE(grid3_ir)
        IF (ALLOCATED(B_grid_ir)) DEALLOCATE(B_grid_ir)
        IF (ALLOCATED(fieldline_endpoints)) DEALLOCATE(fieldline_endpoints)
        IF (ALLOCATED(fieldline_Q)) DEALLOCATE(fieldline_Q)
        IF (ALLOCATED(fieldline_allpos)) DEALLOCATE(fieldline_allpos)
        IF (ALLOCATED(fieldline_connection)) DEALLOCATE(fieldline_connection)
        IF (ALLOCATED(fieldline_user)) DEALLOCATE(fieldline_user)
        IF (ALLOCATED(fieldline_pts)) DEALLOCATE(fieldline_pts)
        IF (ALLOCATED(fieldline_ptn)) DEALLOCATE(fieldline_ptn)

        if (user_defined) then
          call cleanup_user_defined
        end if
      end subroutine cleanup



end module UFiT_Functions_Fortran

