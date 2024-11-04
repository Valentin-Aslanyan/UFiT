

module UFiT_User_Functions
      USE UFiT_Definitions_Fortran
      implicit none


      contains


      subroutine prepare_user_defined
      !Will be executed if running in user mode

        num_ud_variables = 1
        ALLOCATE(fieldline_user(num_ud_variables,numin_tot))

      end subroutine prepare_user_defined


      subroutine trace_cartesian_user(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                      B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                      pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace fieldlines, calculate user defined quantity (Cartesian)

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
        REAL(num) :: B_squared_end1, B_squared_end2

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
        call Bfull_interp(idx_in, pos_next, B_out)
        B_squared_end1 = B_out(1)**2+B_out(2)**2+B_out(3)**2

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
        call Bfull_interp(idx_in, pos_next, B_out)
        B_squared_end2 = B_out(1)**2+B_out(2)**2+B_out(3)**2

        fieldline_user(1,idx_t) = B_squared_end1/B_squared_end2

      end subroutine trace_cartesian_user


      subroutine trace_spherical_user(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                      B_gradB_interp,single_step,pos_start,idx_t,pos_endpoints, &
                                      pos_Q,pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace fieldlines, calculate user defined quantity (spherical)

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
        REAL(num) :: B_squared_end1, B_squared_end2

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
        call Bfull_interp(idx_in, pos_next, B_out)
        B_squared_end1 = B_out(1)**2+B_out(2)**2+B_out(3)**2

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
        call Bfull_interp(idx_in, pos_next, B_out)
        B_squared_end2 = B_out(1)**2+B_out(2)**2+B_out(3)**2

        fieldline_user(1,idx_t) = B_squared_end1/B_squared_end2

      end subroutine trace_spherical_user


      subroutine cleanup_user_defined
      !Will be executed at the end; use this to deallocate, or do any other housekeeping on exit

        IF (ALLOCATED(fieldline_user)) DEALLOCATE(fieldline_user)

      end subroutine cleanup_user_defined


end module UFiT_User_Functions

