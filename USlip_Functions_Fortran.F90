!!
!Coordinate indexing (note: Fortran starts with 1): 
!  1 = X for Cartesian, r for spherical,     R for cylindrical,   log(R) for SphericalExponential/sphex
!  2 = Y for Cartesian, theta for spherical, phi for cylindrical, theta for SphericalExponential/sphex
!  3 = Z for Cartesian, phi for spherical,   Z for cylindrical,   phi for SphericalExponential/sphex
!B(3,grid1,grid2,grid3)
!For irregular (unstructured) grids, there is an addition num_blocks dimension

module USlip_Functions_Fortran
      USE UFiT_Definitions_Fortran
#ifdef _OPENMP
      USE OMP_LIB
#endif
#ifdef USE_NC
      USE netcdf
#endif
      implicit none


      CHARACTER(len=str_mx) :: out_directory
      !Regular coordinate grid (also staggered)
      REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: j_grid
      REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: curlj_grid
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: alpha_grid
      REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: sigma_grid
      REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: sigmaalpha_grid
      !Irregular coordinate grid (e.g. ARMS blocks)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: block_neighbours_ir
      REAL(num), DIMENSION(:,:,:,:,:), ALLOCATABLE :: j_grid_ir
      REAL(num), DIMENSION(:,:,:,:,:), ALLOCATABLE :: curlj_grid_ir
      REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: alpha_grid_ir
      REAL(num), DIMENSION(:,:,:,:,:), ALLOCATABLE :: sigma_grid_ir
      REAL(num), DIMENSION(:,:,:,:,:), ALLOCATABLE :: sigmaalpha_grid_ir

      abstract interface
        subroutine adjacent_points_iface (pt_indices, adj)
          IMPORT
          INTEGER :: pt_indices(4), adj(4,6)
        end subroutine adjacent_points_iface

        subroutine regular_curl_iface (pt_indices, grid_data, curl_out)
          IMPORT
          INTEGER :: pt_indices(3)
          REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: grid_data
          REAL(num) :: curl_out(3)
        end subroutine regular_curl_iface

        subroutine regular_grad_iface (pt_indices, grid_data, grad_out)
          IMPORT
          INTEGER :: pt_indices(3)
          REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: grid_data
          REAL(num) :: grad_out(3)
        end subroutine regular_grad_iface

        subroutine irregular_curl_iface (pt_indices, grid_data, curl_out)
          IMPORT
          INTEGER :: pt_indices(4)
          REAL(num), DIMENSION(:,:,:,:,:), ALLOCATABLE :: grid_data
          REAL(num) :: curl_out(3)
        end subroutine irregular_curl_iface

        subroutine irregular_grad_iface (pt_indices, grid_data, grad_out)
          IMPORT
          INTEGER :: pt_indices(4)
          REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: grid_data
          REAL(num) :: grad_out(3)
        end subroutine irregular_grad_iface
      end interface

      !public :: 


      contains


      subroutine parse_slip_command_args

        INTEGER :: idx, num_args, stat
        CHARACTER(len=str_mx) :: arg

        out_directory = 'out'

        num_args = iargc()

        if (num_args .le. 0) then
          print *, 'You should provide some command line arguments, such as'
          print *, '-b for the input filename'
          call EXIT(136)
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

              case ('-np', '--num_proc', '--NUM_PROC')
                idx = idx + 1
                CALL getarg(idx, arg)
                read(arg,*,iostat=stat) num_proc
                if (stat .ne. 0) then
                  print *, 'unrecognised processor number: ', TRIM(arg)
                end if

              case ('-c', '--command_file')
                idx = idx + 1
                CALL getarg(idx, arg)
                cmd_filename = TRIM(arg)
                read_command_file = .true.

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

              case ('-od', '--output_directory')
                idx = idx + 1
                CALL getarg(idx, arg)
                out_directory = TRIM(arg)

              case default
                print *, 'unrecognised command-line option: ', TRIM(arg)

            end select
            idx = idx + 1
          END DO
        end if

      end subroutine parse_slip_command_args


      subroutine parse_slip_command_file

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

                case ('NP:')
                  read(arg2,*,iostat=stat2) num_proc
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised number of processors: ', TRIM(arg2)
                  end if

                case ('C:')
                  cmd_filename = TRIM(arg2)
                  read_command_file = .true.

                case ('B:')
                  B_filename = TRIM(arg2)

                case ('BT:')
                  read(arg2,*,iostat=stat2) Bfile_type
                  if (stat2 .ne. 0) then
                    print *, 'unrecognised B file type: ', TRIM(arg2)
                  end if

                case ('OD:')
                  out_directory = TRIM(arg2)

                case default
                  print *, 'unrecognised command file option: ', TRIM(arg)

              end select
            ELSE
              continue_read = .false.
            END IF
          END DO

          close(cmd_unit)
        END IF

      end subroutine parse_slip_command_file


      subroutine regular_grad_cartesian(pt_indices, grid_data, grad_out)

        INTEGER :: pt_indices(3)
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: grad_out(3)

        INTEGER :: adj(6)

        call regular_grid_adjacent(pt_indices, periodic_X, periodic_Y, periodic_Z, adj)

        grad_out(1) = (grid_data(adj(1),pt_indices(2),pt_indices(3)) &
                      -grid_data(adj(2),pt_indices(2),pt_indices(3))) &
                      /(grid1(adj(1)) - grid1(adj(2)))

        grad_out(2) = (grid_data(pt_indices(1),adj(3),pt_indices(3)) &
                      -grid_data(pt_indices(1),adj(4),pt_indices(3))) &
                      /(grid2(adj(3)) - grid2(adj(4)))

        grad_out(3) = (grid_data(pt_indices(1),pt_indices(2),adj(5)) &
                      -grid_data(pt_indices(1),pt_indices(2),adj(6))) &
                      /(grid3(adj(5)) - grid3(adj(6)))

      end subroutine regular_grad_cartesian


      subroutine regular_grad_spherical(pt_indices, grid_data, grad_out)

        INTEGER :: pt_indices(3)
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: grad_out(3)

        INTEGER :: adj(6)

        call regular_grid_adjacent(pt_indices, .false., .false., periodic_PHI, adj)

        grad_out(1) = (grid_data(adj(1),pt_indices(2),pt_indices(3)) &
                      -grid_data(adj(2),pt_indices(2),pt_indices(3))) &
                      /(grid1(adj(1)) - grid1(adj(2)))

        grad_out(2) = 2.0_num/(grid1(adj(1)) + grid1(adj(2))) * &
                      (grid_data(pt_indices(1),adj(3),pt_indices(3)) &
                      -grid_data(pt_indices(1),adj(4),pt_indices(3))) &
                      /(grid2(adj(3)) - grid2(adj(4)))

        grad_out(3) = 2.0_num/(grid1(adj(1)) + grid1(adj(2))) / &
                      SIN(0.5_num*(grid2(adj(3)) + grid2(adj(4)))) * &
                      (grid_data(pt_indices(1),pt_indices(2),adj(5)) &
                      -grid_data(pt_indices(1),pt_indices(2),adj(6))) &
                      /(grid3(adj(5)) - grid3(adj(6)))

      end subroutine regular_grad_spherical


      subroutine regular_grad_cylindrical(pt_indices, grid_data, grad_out)

        INTEGER :: pt_indices(3)
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: grad_out(3)

        INTEGER :: adj(6)

        call regular_grid_adjacent(pt_indices, .false., periodic_PHI, periodic_Z, adj)

        grad_out(1) = (grid_data(adj(1),pt_indices(2),pt_indices(3)) &
                      -grid_data(adj(2),pt_indices(2),pt_indices(3))) &
                      /(grid1(adj(1)) - grid1(adj(2)))

        grad_out(2) = 2.0_num/(grid1(adj(1)) + grid1(adj(2))) * &
                      (grid_data(pt_indices(1),adj(3),pt_indices(3)) &
                      -grid_data(pt_indices(1),adj(4),pt_indices(3))) &
                      /(grid2(adj(3)) - grid2(adj(4)))

        grad_out(3) = (grid_data(pt_indices(1),pt_indices(2),adj(5)) &
                      -grid_data(pt_indices(1),pt_indices(2),adj(6))) &
                      /(grid3(adj(5)) - grid3(adj(6)))

      end subroutine regular_grad_cylindrical


      subroutine regular_curl_cartesian(pt_indices, grid_data, curl_out)

        INTEGER :: pt_indices(3)
        REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: curl_out(3)

        INTEGER :: adj(6)

        call regular_grid_adjacent(pt_indices, periodic_X, periodic_Y, periodic_Z, adj)

        curl_out(1) =  (grid_data(3,pt_indices(1),adj(3),pt_indices(3)) &
                       -grid_data(3,pt_indices(1),adj(4),pt_indices(3))) &
                       /(grid2(adj(3)) - grid2(adj(4))) &
                      -(grid_data(2,pt_indices(1),pt_indices(2),adj(5)) &
                       -grid_data(2,pt_indices(1),pt_indices(2),adj(6))) &
                       /(grid3(adj(5)) - grid3(adj(6)))

        curl_out(2) =  (grid_data(1,pt_indices(1),pt_indices(2),adj(5)) &
                       -grid_data(1,pt_indices(1),pt_indices(2),adj(6))) &
                       /(grid3(adj(5)) - grid3(adj(6))) &
                      -(grid_data(3,adj(1),pt_indices(2),pt_indices(3)) &
                       -grid_data(3,adj(2),pt_indices(2),pt_indices(3))) &
                       /(grid1(adj(1)) - grid1(adj(2)))

        curl_out(3) =  (grid_data(2,adj(1),pt_indices(2),pt_indices(3)) &
                       -grid_data(2,adj(2),pt_indices(2),pt_indices(3))) &
                       /(grid1(adj(1)) - grid1(adj(2))) &
                      -(grid_data(1,pt_indices(1),adj(3),pt_indices(3)) &
                       -grid_data(1,pt_indices(1),adj(4),pt_indices(3))) &
                       /(grid2(adj(3)) - grid2(adj(4)))

      end subroutine regular_curl_cartesian


      subroutine regular_curl_spherical(pt_indices, grid_data, curl_out)

        INTEGER :: pt_indices(3)
        REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: curl_out(3)

        INTEGER :: adj(6)

        call regular_grid_adjacent(pt_indices, .false., .false., periodic_PHI, adj)

        curl_out(1) =  2.0_num/(grid1(adj(1)) + grid1(adj(2))) / &
                       SIN(0.5_num*(grid2(adj(3)) + grid2(adj(4)))) * &
                       (grid_data(3,pt_indices(1),adj(3),pt_indices(3))*SIN(grid2(adj(3))) &
                       -grid_data(3,pt_indices(1),adj(4),pt_indices(3))*SIN(grid2(adj(4)))) &
                       /(grid2(adj(3)) - grid2(adj(4))) &
                      -(grid_data(2,pt_indices(1),pt_indices(2),adj(5)) &
                       -grid_data(2,pt_indices(1),pt_indices(2),adj(6))) &
                       /(grid3(adj(5)) - grid3(adj(6)))

        curl_out(2) =  2.0_num/(grid1(adj(1)) + grid1(adj(2))) * &
                       (grid_data(1,pt_indices(1),pt_indices(2),adj(5)) &
                       -grid_data(1,pt_indices(1),pt_indices(2),adj(6))) &
                       /(grid3(adj(5)) - grid3(adj(6)))/SIN(grid2(pt_indices(2))) &
                      -(grid_data(3,adj(1),pt_indices(2),pt_indices(3))*grid1(adj(1)) &
                       -grid_data(3,adj(2),pt_indices(2),pt_indices(3))*grid1(adj(2))) &
                       /(grid1(adj(1)) - grid1(adj(2)))

        curl_out(3) =  2.0_num/(grid1(adj(1)) + grid1(adj(2))) * &
                       (grid_data(2,adj(1),pt_indices(2),pt_indices(3))*grid1(adj(1)) &
                       -grid_data(2,adj(2),pt_indices(2),pt_indices(3))*grid1(adj(2))) &
                       /(grid1(adj(1)) - grid1(adj(2))) &
                      -(grid_data(1,pt_indices(1),adj(3),pt_indices(3)) &
                       -grid_data(1,pt_indices(1),adj(4),pt_indices(3))) &
                       /(grid2(adj(3)) - grid2(adj(4)))

      end subroutine regular_curl_spherical


      subroutine regular_curl_cylindrical(pt_indices, grid_data, curl_out)

        INTEGER :: pt_indices(3)
        REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: curl_out(3)

        INTEGER :: adj(6)

        call regular_grid_adjacent(pt_indices, .false., periodic_PHI, periodic_Z, adj)

        curl_out(1) =  (grid_data(3,pt_indices(1),adj(3),pt_indices(3)) &
                       -grid_data(3,pt_indices(1),adj(4),pt_indices(3))) &
                       /(grid2(adj(3)) - grid2(adj(4)))/grid1(pt_indices(1)) &
                      -(grid_data(2,pt_indices(1),pt_indices(2),adj(5)) &
                       -grid_data(2,pt_indices(1),pt_indices(2),adj(6))) &
                       /(grid3(adj(5)) - grid3(adj(6)))

        curl_out(2) =  (grid_data(1,pt_indices(1),pt_indices(2),adj(5)) &
                       -grid_data(1,pt_indices(1),pt_indices(2),adj(6))) &
                       /(grid3(adj(5)) - grid3(adj(6))) &
                      -(grid_data(3,adj(1),pt_indices(2),pt_indices(3)) &
                       -grid_data(3,adj(2),pt_indices(2),pt_indices(3))) &
                       /(grid1(adj(1)) - grid1(adj(2)))

        curl_out(3) =  2.0_num/(grid1(adj(1)) + grid1(adj(2))) * &
                       (grid_data(2,adj(1),pt_indices(2),pt_indices(3))*grid1(adj(1)) &
                       -grid_data(2,adj(2),pt_indices(2),pt_indices(3))*grid1(adj(2))) &
                       /(grid1(adj(1)) - grid1(adj(2))) &
                      -(grid_data(1,pt_indices(1),adj(3),pt_indices(3)) &
                       -grid_data(1,pt_indices(1),adj(4),pt_indices(3))) &
                       /(grid2(adj(3)) - grid2(adj(4)))

      end subroutine regular_curl_cylindrical


      subroutine regular_grid_adjacent(pt_indices, periodic1, periodic2, periodic3, adj)
      !Each point on the grid is denoted by 3 indices:
      !(1) idx_grid1, (2) idx_grid2, (3) idx_grid3
      !For a given such point with pt_indices, find the 6 adjacent points
      !adj(1,2) has +,- in grid1
      !adj(3,4) has +,- in grid2
      !adj(5,6) has +,- in grid3
      !Return the same point if on edge and not periodic

        INTEGER :: pt_indices(3)
        LOGICAL :: periodic1, periodic2, periodic3
        INTEGER :: adj(6)

        if (pt_indices(1) .eq. 1) then
          if (periodic1) then
            adj(2) = sz_1
          else
            adj(2) = 1
          end if
          adj(1) = pt_indices(1)+1
        else if (pt_indices(1) .eq. sz_1) then
          if (periodic1) then
            adj(1) = 1
          else
            adj(1) = sz_1
          end if
          adj(2) = pt_indices(1)-1
        else
          adj(1) = pt_indices(1)+1
          adj(2) = pt_indices(1)-1
        end if
        if (pt_indices(2) .eq. 1) then
          if (periodic2) then
            adj(4) = sz_2
          else
            adj(4) = 1
          end if
          adj(3) = pt_indices(2)+1
        else if (pt_indices(2) .eq. sz_2) then
          if (periodic2) then
            adj(3) = 1
          else
            adj(3) = sz_2
          end if
          adj(4) = pt_indices(2)-1
        else
          adj(3) = pt_indices(2)+1
          adj(4) = pt_indices(2)-1
        end if
        if (pt_indices(3) .eq. 1) then
          if (periodic3) then
            adj(6) = sz_3
          else
            adj(6) = 1
          end if
          adj(5) = pt_indices(3)+1
        else if (pt_indices(3) .eq. sz_3) then
          if (periodic3) then
            adj(5) = 1
          else
            adj(5) = sz_3
          end if
          adj(6) = pt_indices(3)-1
        else
          adj(5) = pt_indices(3)+1
          adj(6) = pt_indices(3)-1
        end if

      end subroutine regular_grid_adjacent


      subroutine irregular_grad_cartesian(pt_indices, grid_data, grad_out)

        INTEGER :: pt_indices(4)
        REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: grad_out(3)

        INTEGER :: adj(4,6)
        REAL(num) :: c_center(3), c_up(3), c_down(3), c_right(3), c_left(3), c_fwd(3), c_back(3)

        call irregular_grid_adjacent_points(pt_indices, adj)
        call irregular_grid_point_coordinates(pt_indices,c_center)
        call irregular_grid_point_coordinates(adj(:,1),c_up)
        call irregular_grid_point_coordinates(adj(:,2),c_down)
        call irregular_grid_point_coordinates(adj(:,3),c_right)
        call irregular_grid_point_coordinates(adj(:,4),c_left)
        call irregular_grid_point_coordinates(adj(:,5),c_fwd)
        call irregular_grid_point_coordinates(adj(:,6),c_back)

        grad_out(1) = (grid_data(adj(2,1),adj(3,1),adj(4,1),adj(1,1)) &
                      -grid_data(adj(2,2),adj(3,2),adj(4,2),adj(1,2))) &
                      /(c_up(1) - c_down(1))

        grad_out(2) = (grid_data(adj(2,3),adj(3,3),adj(4,3),adj(1,3)) &
                      -grid_data(adj(2,4),adj(3,4),adj(4,4),adj(1,4))) &
                      /(c_right(2) - c_left(2))

        grad_out(3) = (grid_data(adj(2,5),adj(3,5),adj(4,5),adj(1,5)) &
                      -grid_data(adj(2,6),adj(3,6),adj(4,6),adj(1,6))) &
                      /(c_fwd(3) - c_back(3))

      end subroutine irregular_grad_cartesian


      subroutine irregular_grad_spherical(pt_indices, grid_data, grad_out)

        INTEGER :: pt_indices(4)
        REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: grad_out(3)

        INTEGER :: adj(4,6)
        REAL(num) :: c_center(3), c_up(3), c_down(3), c_right(3), c_left(3), c_fwd(3), c_back(3)

        call irregular_grid_adjacent_points(pt_indices, adj)
        call irregular_grid_point_coordinates(pt_indices,c_center)
        call irregular_grid_point_coordinates(adj(:,1),c_up)
        call irregular_grid_point_coordinates(adj(:,2),c_down)
        call irregular_grid_point_coordinates(adj(:,3),c_right)
        call irregular_grid_point_coordinates(adj(:,4),c_left)
        call irregular_grid_point_coordinates(adj(:,5),c_fwd)
        call irregular_grid_point_coordinates(adj(:,6),c_back)

        grad_out(1) = (grid_data(adj(2,1),adj(3,1),adj(4,1),adj(1,1)) &
                      -grid_data(adj(2,2),adj(3,2),adj(4,2),adj(1,2))) &
                      /(c_up(1) - c_down(1))

        grad_out(2) = 2.0_num/(c_up(1) + c_down(1)) * &
                      (grid_data(adj(2,3),adj(3,3),adj(4,3),adj(1,3)) &
                      -grid_data(adj(2,4),adj(3,4),adj(4,4),adj(1,4))) &
                      /(c_right(2) - c_left(2))

        grad_out(3) = 2.0_num/(c_up(1) + c_down(1)) / &
                      SIN(0.5_num*(c_right(2) + c_left(2))) * &
                      (grid_data(adj(2,5),adj(3,5),adj(4,5),adj(1,5)) &
                      -grid_data(adj(2,6),adj(3,6),adj(4,6),adj(1,6))) &
                      /(c_fwd(3) - c_back(3))

      end subroutine irregular_grad_spherical


      subroutine irregular_grad_cylindrical(pt_indices, grid_data, grad_out)

        INTEGER :: pt_indices(4)
        REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: grad_out(3)

        INTEGER :: adj(4,6)
        REAL(num) :: c_center(3), c_up(3), c_down(3), c_right(3), c_left(3), c_fwd(3), c_back(3)

        call irregular_grid_adjacent_points(pt_indices, adj)
        call irregular_grid_point_coordinates(pt_indices,c_center)
        call irregular_grid_point_coordinates(adj(:,1),c_up)
        call irregular_grid_point_coordinates(adj(:,2),c_down)
        call irregular_grid_point_coordinates(adj(:,3),c_right)
        call irregular_grid_point_coordinates(adj(:,4),c_left)
        call irregular_grid_point_coordinates(adj(:,5),c_fwd)
        call irregular_grid_point_coordinates(adj(:,6),c_back)

        grad_out(1) = (grid_data(adj(2,1),adj(3,1),adj(4,1),adj(1,1)) &
                      -grid_data(adj(2,2),adj(3,2),adj(4,2),adj(1,2))) &
                      /(c_up(1) - c_down(1))

        grad_out(2) = 2.0_num/(c_up(1) + c_down(1)) * &
                      (grid_data(adj(2,3),adj(3,3),adj(4,3),adj(1,3)) &
                      -grid_data(adj(2,4),adj(3,4),adj(4,4),adj(1,4))) &
                      /(c_right(2) - c_left(2))

        grad_out(3) = (grid_data(adj(2,5),adj(3,5),adj(4,5),adj(1,5)) &
                      -grid_data(adj(2,6),adj(3,6),adj(4,6),adj(1,6))) &
                      /(c_fwd(3) - c_back(3))

      end subroutine irregular_grad_cylindrical


      subroutine irregular_curl_cartesian(pt_indices, grid_data, curl_out)

        INTEGER :: pt_indices(4)
        REAL(num), DIMENSION(:,:,:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: curl_out(3)

        INTEGER :: adj(4,6)
        REAL(num) :: c_center(3), c_up(3), c_down(3), c_right(3), c_left(3), c_fwd(3), c_back(3)

        call irregular_grid_adjacent_points(pt_indices, adj)
        call irregular_grid_point_coordinates(pt_indices,c_center)
        call irregular_grid_point_coordinates(adj(:,1),c_up)
        call irregular_grid_point_coordinates(adj(:,2),c_down)
        call irregular_grid_point_coordinates(adj(:,3),c_right)
        call irregular_grid_point_coordinates(adj(:,4),c_left)
        call irregular_grid_point_coordinates(adj(:,5),c_fwd)
        call irregular_grid_point_coordinates(adj(:,6),c_back)

        curl_out(1) =  (grid_data(3,adj(2,3),adj(3,3),adj(4,3),adj(1,3)) &
                       -grid_data(3,adj(2,4),adj(3,4),adj(4,4),adj(1,4))) &
                       /(c_right(2) - c_left(2)) &
                      -(grid_data(2,adj(2,5),adj(3,5),adj(4,5),adj(1,5)) &
                       -grid_data(2,adj(2,6),adj(3,6),adj(4,6),adj(1,6))) &
                       /(c_fwd(3) - c_back(3))

        curl_out(2) =  (grid_data(1,adj(2,5),adj(3,5),adj(4,5),adj(1,5)) &
                       -grid_data(1,adj(2,6),adj(3,6),adj(4,6),adj(1,6))) &
                       /(c_fwd(3) - c_back(3)) &
                      -(grid_data(3,adj(2,1),adj(3,1),adj(4,1),adj(1,1)) &
                       -grid_data(3,adj(2,2),adj(3,2),adj(4,2),adj(1,2))) &
                       /(c_up(1) - c_down(1))

        curl_out(3) =  (grid_data(2,adj(2,1),adj(3,1),adj(4,1),adj(1,1)) &
                       -grid_data(2,adj(2,2),adj(3,2),adj(4,2),adj(1,2))) &
                       /(c_up(1) - c_down(1)) &
                      -(grid_data(1,adj(2,3),adj(3,3),adj(4,3),adj(1,3)) &
                       -grid_data(1,adj(2,4),adj(3,4),adj(4,4),adj(1,4))) &
                       /(c_right(2) - c_left(2))

      end subroutine irregular_curl_cartesian


      subroutine irregular_curl_spherical(pt_indices, grid_data, curl_out)

        INTEGER :: pt_indices(4)
        REAL(num), DIMENSION(:,:,:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: curl_out(3)

        INTEGER :: adj(4,6)
        REAL(num) :: c_center(3), c_up(3), c_down(3), c_right(3), c_left(3), c_fwd(3), c_back(3)

        call irregular_grid_adjacent_points(pt_indices, adj)
        call irregular_grid_point_coordinates(pt_indices,c_center)
        call irregular_grid_point_coordinates(adj(:,1),c_up)
        call irregular_grid_point_coordinates(adj(:,2),c_down)
        call irregular_grid_point_coordinates(adj(:,3),c_right)
        call irregular_grid_point_coordinates(adj(:,4),c_left)
        call irregular_grid_point_coordinates(adj(:,5),c_fwd)
        call irregular_grid_point_coordinates(adj(:,6),c_back)

        curl_out(1) =  2.0_num/(c_up(1) + c_down(1)) / &
                       SIN(0.5_num*(c_right(2) + c_left(2))) * &
                      ((grid_data(3,adj(2,3),adj(3,3),adj(4,3),adj(1,3))*SIN(c_right(2)) &
                       -grid_data(3,adj(2,4),adj(3,4),adj(4,4),adj(1,4))*SIN(c_left(2))) &
                       /(c_right(2) - c_left(2)) &
                      -(grid_data(2,adj(2,5),adj(3,5),adj(4,5),adj(1,5)) &
                       -grid_data(2,adj(2,6),adj(3,6),adj(4,6),adj(1,6))) &
                       /(c_fwd(3) - c_back(3)))

        curl_out(2) = 2.0_num/(c_up(1) + c_down(1)) * &
                      ((grid_data(1,adj(2,5),adj(3,5),adj(4,5),adj(1,5)) &
                       -grid_data(1,adj(2,6),adj(3,6),adj(4,6),adj(1,6))) &
                       /(c_fwd(3) - c_back(3))/SIN(c_center(2)) &
                      -(grid_data(3,adj(2,1),adj(3,1),adj(4,1),adj(1,1))*c_up(1) &
                       -grid_data(3,adj(2,2),adj(3,2),adj(4,2),adj(1,2))*c_down(1)) &
                       /(c_up(1) - c_down(1)))

        curl_out(3) = 2.0_num/(c_up(1) + c_down(1)) * &
                      ((grid_data(2,adj(2,1),adj(3,1),adj(4,1),adj(1,1))*c_up(1) &
                       -grid_data(2,adj(2,2),adj(3,2),adj(4,2),adj(1,2))*c_down(1)) &
                       /(c_up(1) - c_down(1)) &
                      -(grid_data(1,adj(2,3),adj(3,3),adj(4,3),adj(1,3)) &
                       -grid_data(1,adj(2,4),adj(3,4),adj(4,4),adj(1,4))) &
                       /(c_right(2) - c_left(2)))

      end subroutine irregular_curl_spherical


      subroutine irregular_curl_cylindrical(pt_indices, grid_data, curl_out)

        INTEGER :: pt_indices(4)
        REAL(num), DIMENSION(:,:,:,:,:), ALLOCATABLE :: grid_data
        REAL(num) :: curl_out(3)

        INTEGER :: adj(4,6)
        REAL(num) :: c_center(3), c_up(3), c_down(3), c_right(3), c_left(3), c_fwd(3), c_back(3)

        call irregular_grid_adjacent_points(pt_indices, adj)
        call irregular_grid_point_coordinates(pt_indices,c_center)
        call irregular_grid_point_coordinates(adj(:,1),c_up)
        call irregular_grid_point_coordinates(adj(:,2),c_down)
        call irregular_grid_point_coordinates(adj(:,3),c_right)
        call irregular_grid_point_coordinates(adj(:,4),c_left)
        call irregular_grid_point_coordinates(adj(:,5),c_fwd)
        call irregular_grid_point_coordinates(adj(:,6),c_back)

        curl_out(1) =  (grid_data(3,adj(2,3),adj(3,3),adj(4,3),adj(1,3)) &
                       -grid_data(3,adj(2,4),adj(3,4),adj(4,4),adj(1,4))) &
                       /(c_right(2) - c_left(2))/c_center(1) &
                      -(grid_data(2,adj(2,5),adj(3,5),adj(4,5),adj(1,5)) &
                       -grid_data(2,adj(2,6),adj(3,6),adj(4,6),adj(1,6))) &
                       /(c_fwd(3) - c_back(3))

        curl_out(2) =  (grid_data(1,adj(2,5),adj(3,5),adj(4,5),adj(1,5)) &
                       -grid_data(1,adj(2,6),adj(3,6),adj(4,6),adj(1,6))) &
                       /(c_fwd(3) - c_back(3)) &
                      -(grid_data(3,adj(2,1),adj(3,1),adj(4,1),adj(1,1)) &
                       -grid_data(3,adj(2,2),adj(3,2),adj(4,2),adj(1,2))) &
                       /(c_up(1) - c_down(1))

        curl_out(3) =  2.0_num/(c_up(1) + c_down(1)) * &
                      ((grid_data(2,adj(2,1),adj(3,1),adj(4,1),adj(1,1))*c_up(1) &
                       -grid_data(2,adj(2,2),adj(3,2),adj(4,2),adj(1,2))*c_down(1)) &
                       /(c_up(1) - c_down(1)) &
                      -(grid_data(1,adj(2,3),adj(3,3),adj(4,3),adj(1,3)) &
                       -grid_data(1,adj(2,4),adj(3,4),adj(4,4),adj(1,4))) &
                       /(c_right(2) - c_left(2)))

      end subroutine irregular_curl_cylindrical


      subroutine compute_block_neighbours(periodic1, periodic2, periodic3)
      !Try to find the indices of neighbouring blocks
      !block_neighbours_ir(1,i) | (2,i) are the +/- blocks in coord1 etc
      !If no such neighbour exists, set the index to -1

        LOGICAL :: periodic1, periodic2, periodic3

        INTEGER :: idx, idx2
        LOGICAL :: neighbour1unset, neighbour2unset, neighbour3unset
        LOGICAL :: neighbour4unset, neighbour5unset, neighbour6unset
        LOGICAL :: coord1aligned, coord2aligned, coord3aligned

        ALLOCATE(block_neighbours_ir(6,num_blocks))
        block_neighbours_ir(:,:) = -1

      !$OMP parallel num_threads(num_proc) private(idx2, neighbour1unset, neighbour2unset, &
      !$OMP& neighbour3unset, neighbour4unset, neighbour5unset, neighbour6unset, coord1aligned, &
      !$OMP& coord2aligned, coord3aligned)
      !$OMP do
        do idx = 1, num_blocks
          neighbour1unset = .true.
          neighbour2unset = .true.
          neighbour3unset = .true.
          neighbour4unset = .true.
          neighbour5unset = .true.
          neighbour6unset = .true.
          idx2 = 1
          do while ((idx2 .le. num_blocks) .and. (neighbour1unset .or. neighbour2unset .or. &
                  neighbour3unset .or. neighbour4unset .or. neighbour5unset .or. neighbour6unset))
            if ((idx .ne. idx2)) then
              coord1aligned = ((grid1_ir(1,idx) .eq. grid1_ir(1,idx2)) .and. (grid1_ir(2,idx) &
                                .eq. grid1_ir(2,idx2)))
              coord2aligned = ((grid2_ir(1,idx) .eq. grid2_ir(1,idx2)) .and. (grid2_ir(2,idx) &
                                .eq. grid2_ir(2,idx2)))
              coord3aligned = ((grid3_ir(1,idx) .eq. grid3_ir(1,idx2)) .and. (grid3_ir(2,idx) &
                                .eq. grid3_ir(2,idx2)))
              if ((grid1_ir(2,idx) .eq. grid1_ir(1,idx2)) .and. coord2aligned .and. &
                   coord3aligned) then
                block_neighbours_ir(1,idx) = idx2
                neighbour1unset = .false.
              else if ((periodic1) .and. (grid1_ir(2,idx) .eq. grid1max) .and. &
                       (grid1_ir(1,idx2) .eq. grid1min)) then
                block_neighbours_ir(1,idx) = idx2
                neighbour1unset = .false.
              end if
              if ((grid1_ir(1,idx) .eq. grid1_ir(2,idx2)) .and. coord2aligned .and. &
                   coord3aligned) then
                block_neighbours_ir(2,idx) = idx2
                neighbour2unset = .false.
              else if ((periodic1) .and. (grid1_ir(1,idx) .eq. grid1min) .and. &
                       (grid1_ir(2,idx2) .eq. grid1max)) then
                block_neighbours_ir(2,idx) = idx2
                neighbour2unset = .false.
              end if
              if ((grid2_ir(2,idx) .eq. grid2_ir(1,idx2)) .and. coord1aligned .and. &
                   coord3aligned) then
                block_neighbours_ir(3,idx) = idx2
                neighbour3unset = .false.
              else if ((periodic2) .and. (grid2_ir(2,idx) .eq. grid2max) .and. &
                       (grid2_ir(1,idx2) .eq. grid2min)) then
                block_neighbours_ir(3,idx) = idx2
                neighbour3unset = .false.
              end if
              if ((grid2_ir(1,idx) .eq. grid2_ir(2,idx2)) .and. coord1aligned .and. &
                   coord3aligned) then
                block_neighbours_ir(4,idx) = idx2
                neighbour4unset = .false.
              else if ((periodic2) .and. (grid2_ir(1,idx) .eq. grid2min) .and. &
                       (grid2_ir(2,idx2) .eq. grid2max)) then
                block_neighbours_ir(4,idx) = idx2
                neighbour4unset = .false.
              end if
              if ((grid3_ir(2,idx) .eq. grid3_ir(1,idx2)) .and. coord1aligned .and. &
                   coord2aligned) then
                block_neighbours_ir(5,idx) = idx2
                neighbour5unset = .false.
              else if ((periodic3) .and. (grid3_ir(2,idx) .eq. grid3max) .and. &
                       (grid3_ir(1,idx2) .eq. grid3min)) then
                block_neighbours_ir(5,idx) = idx2
                neighbour5unset = .false.
              end if
              if ((grid3_ir(1,idx) .eq. grid3_ir(2,idx2)) .and. coord1aligned .and. &
                   coord2aligned) then
                block_neighbours_ir(6,idx) = idx2
                neighbour6unset = .false.
              else if ((periodic3) .and. (grid3_ir(1,idx) .eq. grid3min) .and. &
                       (grid3_ir(2,idx2) .eq. grid3max)) then
                block_neighbours_ir(6,idx) = idx2
                neighbour6unset = .false.
              end if
            end if
            idx2 = idx2+1
          end do
        end do
      !$OMP end do
      !$OMP end parallel

      end subroutine compute_block_neighbours


      subroutine irregular_grid_adjacent_points(pt_indices, adj)
      !Each point on the grid is denoted by 4 indices:
      !(1) idx_block, (2) idx_grid1, (3) idx_grid2, (4) idx_grid3
      !For a given such point with pt_indices, find the 6 adjacent points
      !adj(:,i) has +/- in grid1, +/- in grid2, +/- in grid3
      !Return the same point if no matching point in another block

        INTEGER :: pt_indices(4), adj(4,6)

        if (pt_indices(2) .eq. 1) then
          if (block_neighbours_ir(2,pt_indices(1)) .gt. 0) then
            adj(1,2) = block_neighbours_ir(2,pt_indices(1))
            adj(2,2) = sz_1-1
            adj(3,2) = pt_indices(3)
            adj(4,2) = pt_indices(4)
          else
            adj(1,2) = pt_indices(1)
            adj(2,2) = 1
            adj(3,2) = pt_indices(3)
            adj(4,2) = pt_indices(4)
          end if
          adj(1,1) = pt_indices(1)
          adj(2,1) = 2
          adj(3,1) = pt_indices(3)
          adj(4,1) = pt_indices(4)
        else if (pt_indices(2) .eq. sz_1) then
          if (block_neighbours_ir(1,pt_indices(1)) .gt. 0) then
            adj(1,1) = block_neighbours_ir(1,pt_indices(1))
            adj(2,1) = 2
            adj(3,1) = pt_indices(3)
            adj(4,1) = pt_indices(4)
          else
            adj(1,1) = pt_indices(1)
            adj(2,1) = sz_1
            adj(3,1) = pt_indices(3)
            adj(4,1) = pt_indices(4)
          end if
          adj(1,2) = pt_indices(1)
          adj(2,2) = sz_1 - 1
          adj(3,2) = pt_indices(3)
          adj(4,2) = pt_indices(4)
        else
          adj(1,1) = pt_indices(1)
          adj(2,1) = pt_indices(2) + 1
          adj(3,1) = pt_indices(3)
          adj(4,1) = pt_indices(4)
          adj(1,2) = pt_indices(1)
          adj(2,2) = pt_indices(2) - 1
          adj(3,2) = pt_indices(3)
          adj(4,2) = pt_indices(4)
        end if

        if (pt_indices(3) .eq. 1) then
          if (block_neighbours_ir(4,pt_indices(1)) .gt. 0) then
            adj(1,4) = block_neighbours_ir(4,pt_indices(1))
            adj(2,4) = pt_indices(2)
            adj(3,4) = sz_2-1
            adj(4,4) = pt_indices(4)
          else
            adj(1,4) = pt_indices(1)
            adj(2,4) = pt_indices(2)
            adj(3,4) = 1
            adj(4,4) = pt_indices(4)
          end if
          adj(1,3) = pt_indices(1)
          adj(2,3) = pt_indices(2)
          adj(3,3) = 2
          adj(4,3) = pt_indices(4)
        else if (pt_indices(3) .eq. sz_2) then
          if (block_neighbours_ir(3,pt_indices(1)) .gt. 0) then
            adj(1,3) = block_neighbours_ir(3,pt_indices(1))
            adj(2,3) = pt_indices(2)
            adj(3,3) = 2
            adj(4,3) = pt_indices(4)
          else
            adj(1,3) = pt_indices(1)
            adj(2,3) = pt_indices(2)
            adj(3,3) = sz_2
            adj(4,3) = pt_indices(4)
          end if
          adj(1,4) = pt_indices(1)
          adj(2,4) = pt_indices(2)
          adj(3,4) = sz_2 - 1
          adj(4,4) = pt_indices(4)
        else
          adj(1,3) = pt_indices(1)
          adj(2,3) = pt_indices(2)
          adj(3,3) = pt_indices(3) + 1
          adj(4,3) = pt_indices(4)
          adj(1,4) = pt_indices(1)
          adj(2,4) = pt_indices(2)
          adj(3,4) = pt_indices(3) - 1
          adj(4,4) = pt_indices(4)
        end if

        if (pt_indices(4) .eq. 1) then
          if (block_neighbours_ir(6,pt_indices(1)) .gt. 0) then
            adj(1,6) = block_neighbours_ir(6,pt_indices(1))
            adj(2,6) = pt_indices(2)
            adj(3,6) = pt_indices(3)
            adj(4,6) = sz_3-1
          else
            adj(1,6) = pt_indices(1)
            adj(2,6) = pt_indices(2)
            adj(3,6) = pt_indices(3)
            adj(4,6) = 1
          end if
          adj(1,5) = pt_indices(1)
          adj(2,5) = pt_indices(2)
          adj(3,5) = pt_indices(3)
          adj(4,5) = 2
        else if (pt_indices(4) .eq. sz_3) then
          if (block_neighbours_ir(5,pt_indices(1)) .gt. 0) then
            adj(1,5) = block_neighbours_ir(5,pt_indices(1))
            adj(2,5) = pt_indices(2)
            adj(3,5) = pt_indices(3)
            adj(4,5) = 2
          else
            adj(1,5) = pt_indices(1)
            adj(2,5) = pt_indices(2)
            adj(3,5) = pt_indices(3)
            adj(4,5) = sz_3
          end if
          adj(1,6) = pt_indices(1)
          adj(2,6) = pt_indices(2)
          adj(3,6) = pt_indices(3)
          adj(4,6) = sz_3 - 1
        else
          adj(1,5) = pt_indices(1)
          adj(2,5) = pt_indices(2)
          adj(3,5) = pt_indices(3)
          adj(4,5) = pt_indices(4) + 1
          adj(1,6) = pt_indices(1)
          adj(2,6) = pt_indices(2)
          adj(3,6) = pt_indices(3)
          adj(4,6) = pt_indices(4) - 1
        end if

      end subroutine irregular_grid_adjacent_points


      subroutine irregular_grid_point_coordinates(pt_indices,coords_out)
      !Each point on the grid is denoted by 4 indices:
      !(1) idx_block, (2) idx_grid1, (3) idx_grid2, (4) idx_grid3
      !For a given such point with pt_indices, find
      !the numerical values of those coordinates

        INTEGER :: pt_indices(4)
        REAL(num) :: coords_out(3)

        coords_out(1) = REAL(pt_indices(2),num)/REAL(sz_1-1,num)*(grid1_ir(2,pt_indices(1)) &
                        -grid1_ir(1,pt_indices(1))) + grid1_ir(1,pt_indices(1))
        coords_out(2) = REAL(pt_indices(3),num)/REAL(sz_2-1,num)*(grid2_ir(2,pt_indices(1)) &
                        -grid2_ir(1,pt_indices(1))) + grid2_ir(1,pt_indices(1))
        coords_out(3) = REAL(pt_indices(4),num)/REAL(sz_3-1,num)*(grid3_ir(2,pt_indices(1)) &
                        -grid3_ir(1,pt_indices(1))) + grid3_ir(1,pt_indices(1))

      end subroutine irregular_grid_point_coordinates


      subroutine process_Bfield_slip

        if (grid_regular) then
          IF (grid_separate) THEN
            print *, 'Separate/staggered grid setting not yet implemented for USlip'
            call EXIT(141)
          ELSE !Grid not separate
            ALLOCATE(j_grid(3,sz_1,sz_2,sz_3))
            ALLOCATE(curlj_grid(3,sz_1,sz_2,sz_3))
            ALLOCATE(alpha_grid(sz_1,sz_2,sz_3))
            ALLOCATE(sigma_grid(3,sz_1,sz_2,sz_3))
            ALLOCATE(sigmaalpha_grid(3,sz_1,sz_2,sz_3))
          END IF
        else !Grid irregular
          IF (grid_separate) THEN
            print *, 'Separate/staggered grid setting not yet implemented for USlip'
            call EXIT(142)
          ELSE !Grid not separate
            ALLOCATE(j_grid_ir(3,sz_1,sz_2,sz_3,num_blocks))
            ALLOCATE(curlj_grid_ir(3,sz_1,sz_2,sz_3,num_blocks))
            ALLOCATE(alpha_grid_ir(sz_1,sz_2,sz_3,num_blocks))
            ALLOCATE(sigma_grid_ir(3,sz_1,sz_2,sz_3,num_blocks))
            ALLOCATE(sigmaalpha_grid_ir(3,sz_1,sz_2,sz_3,num_blocks))
            if (geometry .eq. 0) then !Cartesian
              call compute_block_neighbours(periodic_X, periodic_Y, periodic_Z)
            else if (geometry .eq. 1) then !Spherical
              call compute_block_neighbours(.false., .false., periodic_PHI)
            end if
          END IF
        end if

      end subroutine process_Bfield_slip


      subroutine slip_calculation

        INTEGER :: idx, idx1, idx2, idx3, pt_indices3(3), pt_indices4(4)
        REAL(num) :: B_curr(3), mod_B_curr, B_hat_curr(3), curl_curr(3), grad_curr(3), Bdot

        procedure (regular_curl_iface), pointer :: regcrl_ptr => null ()
        procedure (regular_grad_iface), pointer :: reggrd_ptr => null ()
        procedure (irregular_curl_iface), pointer :: ircrl_ptr => null ()
        procedure (irregular_grad_iface), pointer :: irgrd_ptr => null ()

        if (geometry .eq. 0) then !Cartesian
          if (grid_regular) then
            regcrl_ptr => regular_curl_cartesian
            reggrd_ptr => regular_grad_cartesian
          else
            ircrl_ptr => irregular_curl_cartesian
            irgrd_ptr => irregular_grad_cartesian
          end if
        else if (geometry .eq. 1) then !Spherical
          if (grid_regular) then
            regcrl_ptr => regular_curl_spherical
            reggrd_ptr => regular_grad_spherical
          else
            ircrl_ptr => irregular_curl_spherical
            irgrd_ptr => irregular_grad_spherical
          end if
        else if (geometry .eq. 2) then !Cylindrical
          if (grid_regular) then
            regcrl_ptr => regular_curl_cylindrical
            reggrd_ptr => regular_grad_cylindrical
          else
            ircrl_ptr => irregular_curl_cylindrical
            irgrd_ptr => irregular_grad_cylindrical
          end if
        end if


        if (grid_regular) then
          IF (grid_separate) THEN
            print *, 'Separate/staggered grid setting not yet implemented for USlip'
            call EXIT(143)
          ELSE !Grid not separate
            !Pass 1 - get curl(B) and alpha
      !$OMP parallel num_threads(num_proc) private(idx2,idx3,pt_indices3,B_curr,curl_curr)
      !$OMP do
            do idx1 = 1, sz_1
              do idx2 = 1, sz_2
                do idx3 = 1, sz_3
                  pt_indices3(1) = idx1
                  pt_indices3(2) = idx2
                  pt_indices3(3) = idx3
                  B_curr(:) = B_grid(:,idx1,idx2,idx3)
                  call regcrl_ptr(pt_indices3, B_grid, curl_curr)
                  j_grid(:,idx1,idx2,idx3) = curl_curr(:)
                  alpha_grid(idx1,idx2,idx3)=(B_curr(1)*curl_curr(1) &
                         +B_curr(2)*curl_curr(2)+B_curr(3)*curl_curr(3))/(B_curr(1)*B_curr(1) &
                         +B_curr(2)*B_curr(2)+B_curr(3)*B_curr(3))
                end do
              end do
            end do
      !$OMP end do
      !$OMP end parallel

            !Pass 2 - get curl(curl(B)) and grad(alpha)
      !$OMP parallel num_threads(num_proc) private(idx2,idx3,pt_indices3,B_curr,curl_curr, &
      !$OMP& mod_B_curr, grad_curr, B_hat_curr, Bdot)
      !$OMP do
            do idx1 = 1, sz_1
              do idx2 = 1, sz_2
                do idx3 = 1, sz_3
                  pt_indices3(1) = idx1
                  pt_indices3(2) = idx2
                  pt_indices3(3) = idx3
                  B_curr(:) = B_grid(:,idx1,idx2,idx3)
                  mod_B_curr = SQRT(B_curr(1)*B_curr(1)+B_curr(2)*B_curr(2)+B_curr(3)*B_curr(3))
                  B_hat_curr(:) = B_curr(:)/mod_B_curr
                  call regcrl_ptr(pt_indices3, j_grid, curl_curr)
                  call reggrd_ptr(pt_indices3, alpha_grid, grad_curr)
                  Bdot = curl_curr(1)*B_hat_curr(1)+curl_curr(2)*B_hat_curr(2) &
                         +curl_curr(3)*B_hat_curr(3)
                  sigma_grid(:,idx1,idx2,idx3)=-(curl_curr(:)-B_hat_curr(:)*Bdot) &
                                                   /mod_B_curr
                  sigmaalpha_grid(1,idx1,idx2,idx3) = -(grad_curr(2)*B_hat_curr(3) &
                                                      -grad_curr(3)*B_hat_curr(2))
                  sigmaalpha_grid(2,idx1,idx2,idx3) = -(grad_curr(3)*B_hat_curr(1) &
                                                      -grad_curr(1)*B_hat_curr(3))
                  sigmaalpha_grid(3,idx1,idx2,idx3) = -(grad_curr(1)*B_hat_curr(2) &
                                                      -grad_curr(2)*B_hat_curr(1))
                end do
              end do
            end do
      !$OMP end do
      !$OMP end parallel
          END IF
        else !Grid irregular
          IF (grid_separate) THEN
            print *, 'Separate/staggered grid setting not yet implemented for USlip'
            call EXIT(144)
          ELSE !Grid not separate
            !Pass 1 - get curl(B) and alpha
      !$OMP parallel num_threads(num_proc) private(idx1,idx2,idx3,pt_indices4,B_curr,curl_curr)
      !$OMP do
            do idx = 1, num_blocks
              do idx1 = 1, sz_1
                do idx2 = 1, sz_2
                  do idx3 = 1, sz_3
                    pt_indices4(1) = idx
                    pt_indices4(2) = idx1
                    pt_indices4(3) = idx2
                    pt_indices4(4) = idx3
                    B_curr(:) = B_grid_ir(:,idx1,idx2,idx3,idx)
                    call ircrl_ptr(pt_indices4, B_grid_ir, curl_curr)
                    j_grid_ir(:,idx1,idx2,idx3,idx) = curl_curr(:)
                    alpha_grid_ir(idx1,idx2,idx3,idx)=(B_curr(1)*curl_curr(1) &
                           +B_curr(2)*curl_curr(2)+B_curr(3)*curl_curr(3))/(B_curr(1)*B_curr(1) &
                           +B_curr(2)*B_curr(2)+B_curr(3)*B_curr(3))
                  end do
                end do
              end do
            end do
      !$OMP end do
      !$OMP end parallel

            !Pass 2 - get curl(curl(B)) and grad(alpha)
      !$OMP parallel num_threads(num_proc) private(idx1,idx2,idx3,pt_indices4,B_curr,curl_curr, &
      !$OMP& mod_B_curr, grad_curr, B_hat_curr, Bdot)
      !$OMP do
            do idx = 1, num_blocks
              do idx1 = 1, sz_1
                do idx2 = 1, sz_2
                  do idx3 = 1, sz_3
                    pt_indices4(1) = idx
                    pt_indices4(2) = idx1
                    pt_indices4(3) = idx2
                    pt_indices4(4) = idx3
                    B_curr(:) = B_grid_ir(:,idx1,idx2,idx3,idx)
                    mod_B_curr = SQRT(B_curr(1)*B_curr(1)+B_curr(2)*B_curr(2)+B_curr(3)*B_curr(3))
                    B_hat_curr(:) = B_curr(:)/mod_B_curr
                    call ircrl_ptr(pt_indices4, j_grid_ir, curl_curr)
                    call irgrd_ptr(pt_indices4, alpha_grid_ir, grad_curr)
                    Bdot = curl_curr(1)*B_hat_curr(1)+curl_curr(2)*B_hat_curr(2) &
                           +curl_curr(3)*B_hat_curr(3)
                    sigma_grid_ir(:,idx1,idx2,idx3,idx)=-(curl_curr(:)-B_hat_curr(:)*Bdot) &
                                                         /mod_B_curr
                    sigmaalpha_grid_ir(1,idx1,idx2,idx3,idx) = -(grad_curr(2)*B_hat_curr(3) &
                                                            -grad_curr(3)*B_hat_curr(2))
                    sigmaalpha_grid_ir(2,idx1,idx2,idx3,idx) = -(grad_curr(3)*B_hat_curr(1) &
                                                            -grad_curr(1)*B_hat_curr(3))
                    sigmaalpha_grid_ir(3,idx1,idx2,idx3,idx) = -(grad_curr(1)*B_hat_curr(2) &
                                                            -grad_curr(2)*B_hat_curr(1))
                  end do
                end do
              end do
            end do
      !$OMP end do
      !$OMP end parallel
          END IF
        end if

      end subroutine slip_calculation


      subroutine write_slipped_DUMFRIC

#ifdef USE_NC
      !Compiling with NETCDF enabled

        LOGICAL :: stop_found
        INTEGER :: filesize, chunksize, stat_nc, nc_id1, nc_id2, stp_idx
        INTEGER :: dim1id, dim2id, dim3id, realtype, var_id1, var_id2, var_id3
        INTEGER :: var_id4, var_id5, var_id6, var_id7, var_id8, var_id9
        CHARACTER(len=str_mx) :: out_filename
        CHARACTER, DIMENSION(:), ALLOCATABLE :: tempdata

        stop_found = .false.
        stp_idx = str_mx
        DO WHILE ((stp_idx .gt. 0) .and. (.not. stop_found))
          if ((B_filename(stp_idx:stp_idx) .eq. '\') .or. &
              (B_filename(stp_idx:stp_idx) .eq. '/')) then
            stop_found = .true.
          else
            stp_idx = stp_idx - 1
          end if
        END DO
        out_filename=trim(out_directory) // B_filename(stp_idx+1:)

        !Just copy entire file bytewise - ugly!
        open(newunit=nc_id1,file=B_filename,access='stream')
        open(newunit=nc_id2, file=TRIM(out_filename), access = 'stream', status = 'replace')
        chunksize=4096
        ALLOCATE(tempdata(chunksize))
        inquire(file=B_filename, size=filesize)
        DO WHILE (filesize .gt. 0)
          if (chunksize .gt. filesize) then
            chunksize=filesize
            DEALLOCATE(tempdata)
            ALLOCATE(tempdata(chunksize))
          end if
          read(nc_id1) tempdata
          write(nc_id2) tempdata
          filesize = filesize - chunksize
        END DO
        DEALLOCATE(tempdata)
        close(nc_id1)
        close(nc_id2)

        if (num .eq. 8) then
          realtype = NF90_DOUBLE
        else if (num .eq. 4) then
          realtype = NF90_FLOAT
        else
          print *, 'Unrecognized floating point type'
          call EXIT(146)
        end if
        stat_nc = NF90_OPEN(TRIM(out_filename), NF90_WRITE, nc_id2)
        stat_nc = NF90_INQ_VARID(nc_id2, 'r', dim1id)
        stat_nc = NF90_INQ_VARID(nc_id2, 'th', dim2id)
        stat_nc = NF90_INQ_VARID(nc_id2, 'ph', dim3id)
        stat_nc = NF90_REDEF(nc_id2)
        stat_nc = NF90_DEF_VAR(nc_id2, 'CurlBr', realtype, (/ dim1id, dim2id, dim3id /), var_id1)
        stat_nc = NF90_DEF_VAR(nc_id2, 'CurlBt', realtype, (/ dim1id, dim2id, dim3id /), var_id2)
        stat_nc = NF90_DEF_VAR(nc_id2, 'CurlBp', realtype, (/ dim1id, dim2id, dim3id /), var_id3)
        stat_nc = NF90_DEF_VAR(nc_id2, 'Sigmar', realtype, (/ dim1id, dim2id, dim3id /), var_id4)
        stat_nc = NF90_DEF_VAR(nc_id2, 'Sigmat', realtype, (/ dim1id, dim2id, dim3id /), var_id5)
        stat_nc = NF90_DEF_VAR(nc_id2, 'Sigmap', realtype, (/ dim1id, dim2id, dim3id /), var_id6)
        stat_nc = NF90_DEF_VAR(nc_id2, 'SigmaAlphar', realtype, (/ dim1id, dim2id, dim3id /), &
                               var_id7)
        stat_nc = NF90_DEF_VAR(nc_id2, 'SigmaAlphat', realtype, (/ dim1id, dim2id, dim3id /), &
                               var_id8)
        stat_nc = NF90_DEF_VAR(nc_id2, 'SigmaAlphap', realtype, (/ dim1id, dim2id, dim3id /), &
                               var_id9)
        stat_nc = NF90_ENDDEF(nc_id2)
        stat_nc = NF90_PUT_VAR(nc_id2, var_id1, j_grid(1,:,:,:))
        stat_nc = NF90_PUT_VAR(nc_id2, var_id2, j_grid(2,:,:,:))
        stat_nc = NF90_PUT_VAR(nc_id2, var_id3, j_grid(3,:,:,:))
        stat_nc = NF90_PUT_VAR(nc_id2, var_id4, sigma_grid(1,:,:,:))
        stat_nc = NF90_PUT_VAR(nc_id2, var_id5, sigma_grid(2,:,:,:))
        stat_nc = NF90_PUT_VAR(nc_id2, var_id6, sigma_grid(3,:,:,:))
        stat_nc = NF90_PUT_VAR(nc_id2, var_id7, sigmaalpha_grid(1,:,:,:))
        stat_nc = NF90_PUT_VAR(nc_id2, var_id8, sigmaalpha_grid(2,:,:,:))
        stat_nc = NF90_PUT_VAR(nc_id2, var_id9, sigmaalpha_grid(3,:,:,:))
        stat_nc = NF90_CLOSE(nc_id2)

#else
      !Default compilation: NETCDF disabled
        print *, 'Cannot write DUMFRIC file; UFiT was compiled without netCDF'
        print *, 'Ensure netCDF libraries are correct, then run: make USE_NC=True'
        call EXIT(144)
#endif

      end subroutine write_slipped_DUMFRIC


      subroutine write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                        units, d1, d2, d3, data_curr)
      !Block type 3 in Lare3d/CFA SDF

        INTEGER :: L3D_out, d1, d2, d3
        INTEGER(4) :: string_length
        CHARACTER(len=32) :: block_id
        CHARACTER(len=:), ALLOCATABLE :: block_name
        REAL(num), DIMENSION(:,:,:) :: data_curr
        CHARACTER(len=32) :: units

        INTEGER :: block_header_size, idx
        INTEGER(8) :: next_block_location, data_location, data_length
        INTEGER(4) :: blocktype, datatype, ndims, block_info_length, stagger
        INTEGER(4) :: d1_4, d2_4, d3_4
        REAL(8) :: mult
        REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: arr_r8_3d
        CHARACTER(len=32) :: mesh_id

        blocktype = 3
        datatype = 4
        ndims = 3
        block_info_length = 88
        mesh_id = 'grid                            '
        data_length = 8*d1*d2*d3
        d1_4 = d1
        d2_4 = d2
        d3_4 = d3
        mult = 1.0
        stagger = 0
        block_header_size = 160+string_length
        data_location = ftell(L3D_out) + block_header_size
        next_block_location = data_location + data_length

        write(L3D_out) next_block_location
        write(L3D_out) data_location
        write(L3D_out) block_id
        write(L3D_out) data_length
        write(L3D_out) blocktype
        write(L3D_out) datatype
        write(L3D_out) ndims
        write(L3D_out) block_name
        DO idx=len(block_name)+1,string_length
          write(L3D_out) ' '
        END DO
        write(L3D_out) block_info_length
        write(L3D_out) mult
        write(L3D_out) units
        write(L3D_out) mesh_id
        write(L3D_out) d1_4
        write(L3D_out) d2_4
        write(L3D_out) d3_4
        write(L3D_out) stagger

        ALLOCATE(arr_r8_3d(d1,d2,d3))
        arr_r8_3d(:,:,:) = data_curr(:,:,:)
        write(L3D_out) arr_r8_3d
        DEALLOCATE(arr_r8_3d)

      end subroutine write_Lare3d_datablock


      subroutine write_slipped_Lare3d

        LOGICAL :: stop_found
        INTEGER :: L3D_in, L3D_out, infilesize, idx_b, stp_idx, chunksize, bytestoread
        CHARACTER(len=str_mx) :: out_filename
      !Lare3d/Warwick CFSA SDF variables
        INTEGER(1) :: restart_flag, subdomain_file
        INTEGER(4) :: endianness, sdf_version, sdf_revision,summary_size, nblocks
        INTEGER(4) :: block_header_length, step, jobid1, jobid2, string_length, code_io_version
        INTEGER(4) :: geometry_type
        INTEGER(8) :: first_block_location, next_block_location
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
        CHARACTER, DIMENSION(:), ALLOCATABLE :: tempdata
        CHARACTER :: tempbyte
        CHARACTER(len=32) :: units

        stop_found = .false.
        stp_idx = str_mx
        DO WHILE ((stp_idx .gt. 0) .and. (.not. stop_found))
          if ((B_filename(stp_idx:stp_idx) .eq. '\') .or. &
              (B_filename(stp_idx:stp_idx) .eq. '/')) then
            stop_found = .true.
          else
            stp_idx = stp_idx - 1
          end if
        END DO
        out_filename=trim(out_directory) // B_filename(stp_idx+1:)

        inquire(file=B_filename, size=infilesize)
        open(newunit=L3D_in,file=B_filename,access='stream')
        open(newunit=L3D_out,file=out_filename,access='stream')
      !Read header
        read(L3D_in) sdf1_str
        write(L3D_out) sdf1_str
        read(L3D_in) endianness
        write(L3D_out) endianness
        read(L3D_in) sdf_version
        write(L3D_out) sdf_version
        read(L3D_in) sdf_revision
        write(L3D_out) sdf_revision
        read(L3D_in) code_name
        write(L3D_out) code_name
        read(L3D_in) first_block_location
        write(L3D_out) first_block_location
        read(L3D_in) summary_location
        write(L3D_out) summary_location
        read(L3D_in) summary_size
        write(L3D_out) summary_size
        read(L3D_in) nblocks
        write(L3D_out) nblocks+9
        read(L3D_in) block_header_length
        write(L3D_out) block_header_length
        read(L3D_in) step
        write(L3D_out) step
        read(L3D_in) time
        write(L3D_out) time
        read(L3D_in) jobid1
        write(L3D_out) jobid1
        read(L3D_in) jobid2
        write(L3D_out) jobid2
        read(L3D_in) string_length
        write(L3D_out) string_length
        read(L3D_in) code_io_version
        write(L3D_out) code_io_version
        read(L3D_in) restart_flag
        write(L3D_out) restart_flag
        read(L3D_in) subdomain_file
        write(L3D_out) subdomain_file
        ALLOCATE(character(string_length) :: block_name)

      !Copy all existing blocks
        next_block_location = first_block_location
        bytestoread = next_block_location - ftell(L3D_in)
        DO WHILE (bytestoread .gt. 0)
          read(L3D_in) tempbyte
          write(L3D_out) tempbyte
          bytestoread = bytestoread - 1
        END DO
        do idx_b=1,nblocks
          read(L3D_in) next_block_location
          write(L3D_out) next_block_location
          read(L3D_in) data_location
          write(L3D_out) data_location
          read(L3D_in) block_id
          write(L3D_out) block_id
          bytestoread = next_block_location - ftell(L3D_in)
          if (bytestoread .gt. 0) then
            chunksize=MIN(4096,bytestoread)
            ALLOCATE(tempdata(chunksize))
            DO WHILE (bytestoread .gt. 0)
              if (chunksize .gt. bytestoread) then
                chunksize=bytestoread
                DEALLOCATE(tempdata)
                ALLOCATE(tempdata(chunksize))
              end if
              read(L3D_in) tempdata
              write(L3D_out) tempdata
              bytestoread = bytestoread - chunksize
            END DO
            DEALLOCATE(tempdata)
          end if
        end do

        block_id = 'CurlBx                          '
        block_name = 'Magnetic_Field/CurlBx'
        units = 'A/m^3                           '
        call write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                    units, sz_1, sz_2, sz_3, j_grid(1,:,:,:))
        block_id = 'CurlBy                          '
        block_name = 'Magnetic_Field/CurlBy'
        call write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                    units, sz_1, sz_2, sz_3, j_grid(2,:,:,:))
        block_id = 'CurlBz                          '
        block_name = 'Magnetic_Field/CurlBz'
        call write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                    units, sz_1, sz_2, sz_3, j_grid(3,:,:,:))
        block_id = 'Sigmax                          '
        block_name = 'Magnetic_Field/Sigmax'
        units = 'm/s                             '
        call write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                    units, sz_1, sz_2, sz_3, sigma_grid(1,:,:,:))
        block_id = 'Sigmay                          '
        block_name = 'Magnetic_Field/Sigmay'
        call write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                    units, sz_1, sz_2, sz_3, sigma_grid(2,:,:,:))
        block_id = 'Sigmaz                          '
        block_name = 'Magnetic_Field/Sigmaz'
        call write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                    units, sz_1, sz_2, sz_3, sigma_grid(3,:,:,:))
        block_id = 'SigmaAlphax                       '
        block_name = 'Magnetic_Field/SigmaAlphax'
        call write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                    units, sz_1, sz_2, sz_3, sigmaalpha_grid(1,:,:,:))
        block_id = 'SigmaAlphay                       '
        block_name = 'Magnetic_Field/SigmaAlphay'
        call write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                    units, sz_1, sz_2, sz_3, sigmaalpha_grid(2,:,:,:))
        block_id = 'SigmaAlphaz                       '
        block_name = 'Magnetic_Field/SigmaAlphaz'
        call write_Lare3d_datablock(L3D_out, string_length, block_id, block_name, &
                                    units, sz_1, sz_2, sz_3, sigmaalpha_grid(3,:,:,:))

        bytestoread = infilesize - ftell(L3D_in)
        if (bytestoread .gt. 0) then
          chunksize=MIN(4096,bytestoread)
          ALLOCATE(tempdata(chunksize))
          DO WHILE (bytestoread .gt. 0)
            if (chunksize .gt. bytestoread) then
              chunksize=bytestoread
              DEALLOCATE(tempdata)
              ALLOCATE(tempdata(chunksize))
            end if
            read(L3D_in) tempdata
            write(L3D_out) tempdata
            bytestoread = bytestoread - chunksize
          END DO
          DEALLOCATE(tempdata)
        end if

        close(L3D_in)
        close(L3D_out)
        DEALLOCATE(block_name)

      end subroutine write_slipped_Lare3d


      subroutine write_slipped_ARMS_flicks

        INTEGER :: stp_idx, hdrin_unit, flicksin_unit, hdrout_unit, flicksout_unit, stat, stat2
        INTEGER :: ftrin_unit, ftrout_unit, qtynum
        INTEGER :: idx_leaf, idx_blk, idx1, idx2, idx3, real_size
        LOGICAL :: log_grid, stop_found, bfile_exists
        INTEGER(4) :: ntblks, nlblks, newgrd, nvar
        INTEGER(4) :: iputwrk(21)
        REAL(4) :: time32
        REAL(4) :: rputwrk32(6)
        REAL(4), DIMENSION(:), ALLOCATABLE :: datain_temp32
        REAL(4), DIMENSION(:), ALLOCATABLE :: dataout_temp32
        REAL(8) :: time64
        REAL(8) :: rputwrk64(6)
        REAL(8), DIMENSION(:), ALLOCATABLE :: datain_temp64
        REAL(8), DIMENSION(:), ALLOCATABLE :: dataout_temp64
        REAL(num) :: minvar, maxvar
        CHARACTER(len=str_mx) :: hdr_line
        CHARACTER(len=str_mx) :: ftr_line
        CHARACTER(len=4) :: blank4
        CHARACTER(len=8) :: blank8
        CHARACTER(len=25) :: blank25
        CHARACTER(len=str_mx) :: hdrin_filename
        CHARACTER(len=str_mx) :: ftrin_filename
        CHARACTER(len=str_mx) :: hdrout_filename
        CHARACTER(len=str_mx) :: ftrout_filename
        CHARACTER(len=str_mx) :: out_filename

        stop_found = .false.
        nvar = 0

        stp_idx = str_mx
        DO WHILE ((stp_idx .gt. 0) .and. (.not. stop_found))
          if (B_filename(stp_idx:stp_idx) .eq. '.') then
            stop_found = .true.
          else
            stp_idx = stp_idx - 1
          end if
        END DO
        hdrin_filename=B_filename(1:stp_idx) // 'hdr'
        ftrin_filename=B_filename(1:stp_idx) // 'ftr'
        hdrout_filename=trim(out_directory) // 'flicks.hdr'
        ftrout_filename=trim(out_directory) // 'flicks.ftr'
        out_filename=trim(out_directory) // 'flicks' // B_filename(stp_idx:)

        inquire(file=TRIM(hdrin_filename), exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open flicks header file (must end with .hdr)'
          print *, 'This should be in the same directory as the target flicks file'
          print *, 'Attemped to read filename: '
          print *, TRIM(hdrin_filename)
          call EXIT(130)
        END IF
        open(newunit=hdrin_unit,file=TRIM(hdrin_filename),form="formatted")
        open(newunit=hdrout_unit,file=TRIM(hdrout_filename),form="formatted")
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat) trim(hdr_line)
        if (hdr_line(6:7) .eq. '32') then
          real_size = 4
        else if (hdr_line(6:7) .eq. '64') then
          real_size = 8
        else
          print *, 'Following line could not be interpreted as floating point size:'
          print *, hdr_line
          call EXIT(134)
        end if
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        do while (stat .eq. 0)
          write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
          read(hdr_line(1:1),*,iostat=stat2) qtynum
          nvar = nvar + qtynum
          read(hdrin_unit, '(A)', iostat=stat) hdr_line
        end do
        close(hdrin_unit)
        write(hdrout_unit, '(A)', iostat=stat2) '3CurlB'
        write(hdrout_unit, '(A)', iostat=stat2) '3Sigma'
        write(hdrout_unit, '(A)', iostat=stat2) '3SigmaAlpha'
        close(hdrout_unit)

        inquire(file=TRIM(ftrin_filename), exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open flicks header file (must end with .hdr)'
          print *, 'This should be in the same directory as the target flicks file'
          print *, 'Attemped to read filename: '
          print *, TRIM(hdrin_filename)
          call EXIT(130)
        END IF
        open(newunit=ftrin_unit,file=TRIM(ftrin_filename),form="formatted")
        open(newunit=ftrout_unit,file=TRIM(ftrout_filename),form="formatted")
        read(ftrin_unit, '(A)', iostat=stat) ftr_line
        do while (stat .eq. 0)
          write(ftrout_unit, '(A)', iostat=stat2) trim(ftr_line)
          read(ftrin_unit, '(A)', iostat=stat) ftr_line
        end do
        close(ftrin_unit)
        minvar=MINVAL(j_grid_ir(1,:,:,:,:))
        maxvar=MAXVAL(j_grid_ir(1,:,:,:,:))
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) minvar, 'Min R CurlB'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) maxvar, 'Max R CurlB'
        minvar=MINVAL(j_grid_ir(2,:,:,:,:))
        maxvar=MAXVAL(j_grid_ir(2,:,:,:,:))
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) minvar, 'Min T CurlB'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) maxvar, 'Max T CurlB'
        minvar=MINVAL(j_grid_ir(3,:,:,:,:))
        maxvar=MAXVAL(j_grid_ir(3,:,:,:,:))
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) minvar, 'Min P CurlB'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) maxvar, 'Max P CurlB'
        minvar=MINVAL(j_grid_ir(1,:,:,:,:)**2+j_grid_ir(2,:,:,:,:)**2+j_grid_ir(3,:,:,:,:)**2)
        maxvar=MAXVAL(j_grid_ir(1,:,:,:,:)**2+j_grid_ir(2,:,:,:,:)**2+j_grid_ir(3,:,:,:,:)**2)
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) SQRT(minvar), 'Min CurlB Magnitude'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) SQRT(maxvar), 'Max CurlB Magnitude'
        minvar=MINVAL(sigma_grid_ir(1,:,:,:,:))
        maxvar=MAXVAL(sigma_grid_ir(1,:,:,:,:))
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) minvar, 'Min R Sigma'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) maxvar, 'Max R Sigma'
        minvar=MINVAL(sigma_grid_ir(2,:,:,:,:))
        maxvar=MAXVAL(sigma_grid_ir(2,:,:,:,:))
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) minvar, 'Min T Sigma'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) maxvar, 'Max T Sigma'
        minvar=MINVAL(sigma_grid_ir(3,:,:,:,:))
        maxvar=MAXVAL(sigma_grid_ir(3,:,:,:,:))
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) minvar, 'Min P Sigma'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) maxvar, 'Max P Sigma'
        minvar=MINVAL(sigma_grid_ir(1,:,:,:,:)**2+sigma_grid_ir(2,:,:,:,:)**2 &
                     +sigma_grid_ir(3,:,:,:,:)**2)
        maxvar=MAXVAL(sigma_grid_ir(1,:,:,:,:)**2+sigma_grid_ir(2,:,:,:,:)**2 &
                     +sigma_grid_ir(3,:,:,:,:)**2)
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) SQRT(minvar), 'Min Sigma Magnitude'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) SQRT(maxvar), 'Max Sigma Magnitude'
        minvar=MINVAL(sigmaalpha_grid_ir(1,:,:,:,:))
        maxvar=MAXVAL(sigmaalpha_grid_ir(1,:,:,:,:))
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) minvar, 'Min R SigmaAlpha'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) maxvar, 'Max R SigmaAlpha'
        minvar=MINVAL(sigmaalpha_grid_ir(2,:,:,:,:))
        maxvar=MAXVAL(sigmaalpha_grid_ir(2,:,:,:,:))
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) minvar, 'Min T SigmaAlpha'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) maxvar, 'Max T SigmaAlpha'
        minvar=MINVAL(sigmaalpha_grid_ir(3,:,:,:,:))
        maxvar=MAXVAL(sigmaalpha_grid_ir(3,:,:,:,:))
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) minvar, 'Min P SigmaAlpha'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) maxvar, 'Max P SigmaAlpha'
        minvar=MINVAL(sigmaalpha_grid_ir(1,:,:,:,:)**2+sigmaalpha_grid_ir(2,:,:,:,:)**2 &
                     +sigmaalpha_grid_ir(3,:,:,:,:)**2)
        maxvar=MAXVAL(sigmaalpha_grid_ir(1,:,:,:,:)**2+sigmaalpha_grid_ir(2,:,:,:,:)**2 &
                     +sigmaalpha_grid_ir(3,:,:,:,:)**2)
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) SQRT(minvar), 'Min SigmaAlpha Magnitude'
        write(ftrout_unit, '(E11.4,A)', iostat=stat2) SQRT(maxvar), 'Max SigmaAlpha Magnitude'
        close(ftrout_unit)

        if (real_size .eq. 4) then
          ALLOCATE(datain_temp32(nvar))
          ALLOCATE(dataout_temp32(nvar+9))
        else if (real_size .eq. 8) then
          ALLOCATE(dataout_temp64(nvar))
          ALLOCATE(dataout_temp64(nvar+9))
        end if

        inquire(file=B_filename, exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open flicks file'
          print *, 'Attemped to read filename: '
          print *, TRIM(hdrin_filename)
          call EXIT(133)
        END IF
        open(newunit=flicksin_unit,file=B_filename,access='stream',convert='big_endian')
        open(newunit=flicksout_unit,file=out_filename,access='stream',convert='big_endian')
        read(flicksin_unit) blank25
        write(flicksout_unit) blank25
        if (real_size .eq. 4) then
          read(flicksin_unit) time32
          write(flicksout_unit) time32
        else if (real_size .eq. 8) then
          read(flicksin_unit) time64
          write(flicksout_unit) time64
        end if
        read(flicksin_unit) blank8
        write(flicksout_unit) blank8
        read(flicksin_unit) ntblks
        write(flicksout_unit) ntblks
        read(flicksin_unit) nlblks
        write(flicksout_unit) nlblks
        num_blocks = nlblks
        read(flicksin_unit) newgrd
        write(flicksout_unit) newgrd
        read(flicksin_unit) blank4
        write(flicksout_unit) blank4

        if (real_size .eq. 4) then
          idx_leaf = 1
          DO idx_blk = 1,ntblks
            read(flicksin_unit) blank4
            write(flicksout_unit) blank4
            read(flicksin_unit) iputwrk
            write(flicksout_unit) iputwrk
            read(flicksin_unit) blank8
            write(flicksout_unit) blank8
            read(flicksin_unit) rputwrk32
            write(flicksout_unit) rputwrk32
            read(flicksin_unit) blank4
            write(flicksout_unit) blank4
            if (iputwrk(3) .eq. 1) then
              DO idx3 = 1,sz_3
                DO idx2 = 1,sz_2
                  DO idx1 = 1,sz_1
                    read(flicksin_unit) qtynum
                    write(flicksout_unit) qtynum+36
                    read(flicksin_unit) datain_temp32
                    dataout_temp32(1:nvar)=datain_temp32(:)
                    dataout_temp32(nvar+1:nvar+3)=j_grid_ir(:,idx1,idx2,idx3,idx_leaf)
                    dataout_temp32(nvar+4:nvar+6)=sigma_grid_ir(:,idx1,idx2,idx3,idx_leaf)
                    dataout_temp32(nvar+7:nvar+9)=sigmaalpha_grid_ir(:,idx1,idx2,idx3,idx_leaf)
                    write(flicksout_unit) dataout_temp32
                    read(flicksin_unit) qtynum
                    write(flicksout_unit) qtynum+36
                  END DO
                END DO
              END DO
              idx_leaf = idx_leaf + 1
            end if
          END DO
        else if (real_size .eq. 8) then
          idx_leaf = 1
          DO idx_blk = 1,ntblks
            read(flicksin_unit) blank4
            write(flicksout_unit) blank4
            read(flicksin_unit) iputwrk
            write(flicksout_unit) iputwrk
            read(flicksin_unit) blank8
            write(flicksout_unit) blank8
            read(flicksin_unit) rputwrk64
            write(flicksout_unit) rputwrk64
            read(flicksin_unit) blank4
            write(flicksout_unit) blank4
            if (iputwrk(3) .eq. 1) then
              DO idx3 = 1,sz_3
                DO idx2 = 1,sz_2
                  DO idx1 = 1,sz_1
                    read(flicksin_unit) qtynum
                    write(flicksout_unit) qtynum+72
                    read(flicksin_unit) datain_temp64
                    dataout_temp64(1:nvar)=datain_temp64(:)
                    dataout_temp64(nvar+1:nvar+3)=j_grid_ir(:,idx1,idx2,idx3,idx_leaf)
                    dataout_temp64(nvar+4:nvar+6)=sigma_grid_ir(:,idx1,idx2,idx3,idx_leaf)
                    dataout_temp64(nvar+7:nvar+9)=sigmaalpha_grid_ir(:,idx1,idx2,idx3,idx_leaf)
                    write(flicksout_unit) dataout_temp64
                    read(flicksin_unit) qtynum
                    write(flicksout_unit) qtynum+72
                  END DO
                END DO
              END DO
              idx_leaf = idx_leaf + 1
            end if
          END DO
        end if

        close(flicksin_unit)
        close(flicksout_unit)
        if (real_size .eq. 4) then
          DEALLOCATE(datain_temp32)
          DEALLOCATE(dataout_temp32)
        else if (real_size .eq. 8) then
          DEALLOCATE(datain_temp64)
          DEALLOCATE(dataout_temp64)
        end if

      end subroutine write_slipped_ARMS_flicks


      subroutine write_slipped_ARMS_bfield

        INTEGER :: stp_idx, hdrin_unit, bfieldin_unit, hdrout_unit, bfieldout_unit, stat, stat2
        INTEGER :: qtynum, idx_leaf, idx_blk, idx1, idx2, idx3
        LOGICAL :: log_grid, stop_found, bfile_exists
        INTEGER(4) :: ntblks, nlblks
        INTEGER(4) :: iputwrk(35)
        REAL(8) :: time64
        REAL(8) :: rputwrk64(6)
        REAL(8), DIMENSION(:), ALLOCATABLE :: datain_temp64
        REAL(8), DIMENSION(:), ALLOCATABLE :: dataout_temp64
        CHARACTER(len=str_mx) :: hdr_line
        CHARACTER(len=4) :: blank4
        CHARACTER(len=8) :: blank8
        CHARACTER(len=str_mx) :: hdrin_filename
        CHARACTER(len=str_mx) :: hdrout_filename
        CHARACTER(len=str_mx) :: out_filename

        stop_found = .false.

        stp_idx = str_mx
        DO WHILE ((stp_idx .gt. 0) .and. (.not. stop_found))
          if (B_filename(stp_idx:stp_idx) .eq. '.') then
            stop_found = .true.
          else
            stp_idx = stp_idx - 1
          end if
        END DO
        hdrin_filename=B_filename(1:stp_idx) // 'hdr'
        hdrout_filename=trim(out_directory) // 'bfield.hdr'
        out_filename=trim(out_directory) // 'bfield' // B_filename(stp_idx:)

        inquire(file=TRIM(hdrin_filename), exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open bfield header file (must end with .hdr)'
          print *, 'This should be in the same directory as the target bfield file'
          print *, 'Attemped to read filename: '
          print *, TRIM(hdrin_filename)
          call EXIT(135)
        END IF
        open(newunit=hdrin_unit,file=TRIM(hdrin_filename),form="formatted")
        open(newunit=hdrout_unit,file=TRIM(hdrout_filename),form="formatted")
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        read(hdrin_unit, '(A)', iostat=stat) hdr_line
        write(hdrout_unit, '(A)', iostat=stat2) trim(hdr_line)
        close(hdrin_unit)
        close(hdrout_unit)

        ALLOCATE(datain_temp64(3))
        ALLOCATE(dataout_temp64(12))
        inquire(file=B_filename, exist=bfile_exists)
        IF (.not. bfile_exists) THEN
          print *, 'Unable to open bfield file'
          print *, 'Attemped to read filename: '
          print *, TRIM(B_filename)
          call EXIT(137)
        END IF
        open(newunit=bfieldin_unit,file=B_filename,access='stream',convert='big_endian')
        open(newunit=bfieldout_unit,file=out_filename,access='stream',convert='big_endian')
        read(bfieldin_unit) blank4
        write(bfieldout_unit) blank4
        read(bfieldin_unit) time64
        write(bfieldout_unit) time64
        read(bfieldin_unit) blank8
        write(bfieldout_unit) blank8
        read(bfieldin_unit) ntblks
        write(bfieldout_unit) ntblks
        read(bfieldin_unit) nlblks
        write(bfieldout_unit) nlblks
        num_blocks = nlblks
        read(bfieldin_unit) blank4
        write(bfieldout_unit) blank4

        idx_leaf = 1
        DO idx_blk = 1,ntblks
          read(bfieldin_unit) blank4
          write(bfieldout_unit) blank4
          read(bfieldin_unit) iputwrk
          write(bfieldout_unit) iputwrk
          read(bfieldin_unit) blank8
          write(bfieldout_unit) blank8
          read(bfieldin_unit) rputwrk64
          write(bfieldout_unit) rputwrk64
          read(bfieldin_unit) blank4
          write(bfieldout_unit) blank4
          if (iputwrk(3) .eq. 1) then
            DO idx3 = 1,sz_3
              DO idx2 = 1,sz_2
                DO idx1 = 1,sz_1
                  read(bfieldin_unit) qtynum
                  write(bfieldin_unit) qtynum+72
                  read(bfieldin_unit) datain_temp64
                  dataout_temp64(1:3)=datain_temp64(:)
                  dataout_temp64(4:6)=j_grid_ir(:,idx1,idx2,idx3,idx_leaf)
                  dataout_temp64(7:9)=sigma_grid_ir(:,idx1,idx2,idx3,idx_leaf)
                  dataout_temp64(10:12)=sigmaalpha_grid_ir(:,idx1,idx2,idx3,idx_leaf)
                  read(bfieldin_unit) qtynum
                  write(bfieldin_unit) qtynum+72
                END DO
              END DO
            END DO
            idx_leaf = idx_leaf + 1
          end if
        END DO

        close(bfieldin_unit)
        close(bfieldout_unit)
        DEALLOCATE(datain_temp64)
        DEALLOCATE(dataout_temp64)

      end subroutine write_slipped_ARMS_bfield


      subroutine write_slip_output

        LOGICAL :: stop_found
        INTEGER :: out_unit, stp_idx
        CHARACTER(len=str_mx) :: out_filename

        IF (Bfile_type_actual .eq. 0) THEN
          stop_found = .false.
          stp_idx = str_mx
          DO WHILE ((stp_idx .gt. 0) .and. (.not. stop_found))
            if ((B_filename(stp_idx:stp_idx) .eq. '\') .or. &
                (B_filename(stp_idx:stp_idx) .eq. '/')) then
              stop_found = .true.
            else
              stp_idx = stp_idx - 1
            end if
          END DO
          out_filename=trim(out_directory) // B_filename(stp_idx+1:)

          open(newunit=out_unit,file=out_filename,access='stream')
          write(out_unit) sz_1
          write(out_unit) sz_2
          write(out_unit) sz_3
          write(out_unit) grid1
          write(out_unit) grid2
          write(out_unit) grid3
          write(out_unit) B_grid
          write(out_unit) j_grid
          write(out_unit) sigma_grid
          write(out_unit) sigmaalpha_grid
          close(out_unit)
        ELSE IF (Bfile_type_actual .eq. 10) THEN
          call write_slipped_DUMFRIC
        ELSE IF (Bfile_type_actual .eq. 20) THEN
          call write_slipped_Lare3d
        ELSE IF (Bfile_type_actual .eq. 30) THEN
          call write_slipped_ARMS_flicks
        ELSE IF (Bfile_type_actual .eq. 31) THEN
          call write_slipped_ARMS_bfield
        END IF

      end subroutine write_slip_output


      subroutine slip_cleanup

        IF (ALLOCATED(j_grid)) DEALLOCATE(j_grid)
        IF (ALLOCATED(curlj_grid)) DEALLOCATE(curlj_grid)
        IF (ALLOCATED(alpha_grid)) DEALLOCATE(alpha_grid)
        IF (ALLOCATED(sigma_grid)) DEALLOCATE(sigma_grid)
        IF (ALLOCATED(sigmaalpha_grid)) DEALLOCATE(sigmaalpha_grid)
        IF (ALLOCATED(block_neighbours_ir)) DEALLOCATE(block_neighbours_ir)
        IF (ALLOCATED(j_grid_ir)) DEALLOCATE(j_grid_ir)
        IF (ALLOCATED(curlj_grid_ir)) DEALLOCATE(curlj_grid_ir)
        IF (ALLOCATED(alpha_grid_ir)) DEALLOCATE(alpha_grid_ir)
        IF (ALLOCATED(sigma_grid_ir)) DEALLOCATE(sigma_grid_ir)
        IF (ALLOCATED(sigmaalpha_grid_ir)) DEALLOCATE(sigmaalpha_grid_ir)

      end subroutine slip_cleanup


end module USlip_Functions_Fortran

