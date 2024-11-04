
module UFiT_Definitions_Fortran
      implicit none


      INTEGER, PARAMETER :: num = KIND(0.d0) ! precision for floats
      INTEGER, PARAMETER :: HEADER_LENGTH = 1024
      REAL(num), PARAMETER :: PI = 4.0_num*ATAN(1.0_num)
      REAL(num), PARAMETER :: TWOPI = 8.0_num*ATAN(1.0_num)
      REAL(num), PARAMETER :: EPS = 1.0E-15_num  !Small value
      REAL(num), PARAMETER :: theta_periodicity_threshold = 0.03_num  !If theta grid is within this distance of pole, assume it's open
      INTEGER, PARAMETER :: str_mx = 300     ! max length of character strings
      INTEGER :: geometry, Bfile_type, input_type, numin_tot, numin1, numin2, numin3
      LOGICAL :: save_endpoints, save_Q, save_fieldlines, save_connection, user_defined, check_starts
      LOGICAL :: normalized_B, include_curvature, periodic_X, periodic_Y, periodic_Z, periodic_PHI
      INTEGER :: num_proc
      INTEGER :: MAX_STEPS
      REAL(num) :: step_size
      LOGICAL :: read_command_file, print_devices
      INTEGER :: sz_1, sz_2, sz_3, num_blocks
      CHARACTER(len=str_mx) :: cmd_filename
      CHARACTER(len=str_mx) :: B_filename
      CHARACTER(len=str_mx) :: in_filename
      CHARACTER(len=str_mx) :: out_filename
      LOGICAL :: grid_regular, grid_separate
      !Following are seeds of fieldlines
      REAL(num), DIMENSION(:), ALLOCATABLE :: coord1_in  !X for Cartesian, r for spherical
      REAL(num), DIMENSION(:), ALLOCATABLE :: coord2_in  !Y for Cartesian, theta for spherical
      REAL(num), DIMENSION(:), ALLOCATABLE :: coord3_in  !Z for Cartesian, phi for spherical
      REAL(num) :: grid1min, grid1max, grid2min, grid2max, grid3min, grid3max
      !Regular coordinate grid
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid1      !X for Cartesian, r for spherical
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid2      !Y for Cartesian, theta for spherical
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid3      !Z for Cartesian, phi for spherical
      REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: B_grid
      !Regular coordinate grids, staggered (independent grids for each B component)
      INTEGER :: sz_11, sz_12, sz_13, sz_21, sz_22, sz_23, sz_31, sz_32, sz_33
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid1_1      !X/r coordinates for B_X/B_r
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid1_2      !Y/theta coordinates for B_X/B_r
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid1_3      !Z/phi coordinates for B_X/B_r
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid2_1      !X/r coordinates for B_Y/B_theta
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid2_2      !Y/theta coordinates for B_Y/B_theta
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid2_3      !Z/phi coordinates for B_Y/B_theta
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid3_1      !X/r coordinates for B_Z/B_phi
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid3_2      !Y/theta coordinates for B_Z/B_phi
      REAL(num), DIMENSION(:), ALLOCATABLE :: grid3_3      !Z/phi coordinates for B_Z/B_phi
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: B_grid1  !B_X/B_r
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: B_grid2  !B_Y/B_theta
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: B_grid3  !B_Z/B_phi
      !Irregular coordinate grid (e.g. ARMS blocks)
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: grid1_ir      !X for Cartesian, r for spherical
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: grid2_ir      !Y for Cartesian, theta for spherical
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: grid3_ir      !Z for Cartesian, phi for spherical
      REAL(num), DIMENSION(:,:,:,:,:), ALLOCATABLE :: B_grid_ir
      REAL(num) :: coord_width(3)                        !Distance between last and first points in respective coordinate
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: fieldline_endpoints
      REAL(num), DIMENSION(:), ALLOCATABLE :: fieldline_Q
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: fieldline_allpos
      BYTE, DIMENSION(:), ALLOCATABLE :: fieldline_connection
      REAL(num), DIMENSION(:,:), ALLOCATABLE :: fieldline_user
      INTEGER, DIMENSION(:), ALLOCATABLE :: fieldline_pts
      INTEGER, DIMENSION(:), ALLOCATABLE :: fieldline_ptn
      !Fieldlines are considered closed if both endpoints are within this constant *size(Z) of Z_min or *size(R) of R_min respectively
      REAL(num) :: closed_fl_constant=0.1_num
      REAL(num) :: closed_fl_size
      INTEGER :: num_ud_variables

      abstract interface
        subroutine get_position0_iface (pos_start,idx_in)
          IMPORT
          REAL(num) :: pos_start(3)
          INTEGER :: idx_in
        end subroutine get_position0_iface

        subroutine check_position_iface (pos_in,pos_out,res)
          IMPORT
          REAL(num) :: pos_in(3), pos_out(3)
          LOGICAL :: res
        end subroutine check_position_iface

        subroutine intercept_boundary_iface (pos_in,direction,delta_out)
          IMPORT
          REAL(num) :: pos_in(3), direction(3), delta_out
        end subroutine intercept_boundary_iface

        subroutine B_interp_iface (idx_in, pos_in, B_out)
          IMPORT
          INTEGER :: idx_in(9)
          REAL(num), INTENT(IN) :: pos_in(3)
          REAL(num), INTENT(OUT) :: B_out(3)
        end subroutine B_interp_iface

        subroutine Bfull_interp_iface (idx_in, pos_in, B_out)
          IMPORT
          INTEGER :: idx_in(9)
          REAL(num), INTENT(IN) :: pos_in(3)
          REAL(num), INTENT(OUT) :: B_out(3)
        end subroutine Bfull_interp_iface

        subroutine B_gradB_interp_iface (idx_in, pos_in, B_out, gradB_out)
          IMPORT
          INTEGER :: idx_in(9)
          REAL(num), INTENT(IN) :: pos_in(3)
          REAL(num), INTENT(OUT) :: B_out(3)
          REAL(num), INTENT(OUT) :: gradB_out(3,3)
        end subroutine B_gradB_interp_iface

        subroutine single_step_iface (idx_in, pos_in, pos_out, u_vec, v_vec, dl, keep_running, &
                                  check_position, B_interp, B_gradB_interp)
          IMPORT
          INTEGER :: idx_in(9)
          REAL(num) :: pos_in(3), pos_out(3), u_vec(3), v_vec(3), dl
          LOGICAL :: keep_running
          procedure(check_position_iface) :: check_position
          procedure(B_interp_iface) :: B_interp
          procedure(B_gradB_interp_iface) :: B_gradB_interp
        end subroutine single_step_iface

        subroutine trace_fl_iface (check_position,intercept_boundary,B_interp,Bfull_interp, &
                                B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
          IMPORT
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
        end subroutine trace_fl_iface
      end interface


end module UFiT_Definitions_Fortran

