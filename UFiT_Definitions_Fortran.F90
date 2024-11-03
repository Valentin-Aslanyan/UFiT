
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


end module UFiT_Definitions_Fortran

