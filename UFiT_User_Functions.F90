

!Note that the three trace_XXX_user functions are identical in the GitHub repo,
!but these are placeholders to allow the user to customize different coordinate
!systems separately
module UFiT_User_Functions
      USE UFiT_Definitions_Fortran
      implicit none


      contains


      subroutine prepare_user_defined
      !Will be executed if running in user mode

        num_ud_variables = 1
        ALLOCATE(fieldline_user(num_ud_variables,numin_tot))

      end subroutine prepare_user_defined


      subroutine fieldline_classification_user(idx_t,pos_endpoints,connection_type,Q_sign)
      !This is switched on automatically by -ud/--user_defined flag,
      !OR by -ue/--user_endpoints if you want to use only this routine, but trace normally otherwise

      !Use this to define the "type" of endpoint, e.g. open/closed/disconnected (can be integer [-128,127])
      !Also a multiplier for Q (suggested: +/-1.0_num for closed/open)
      !pos_endpoints(:,idx_t) contains the (coord1, coord2, coord3) location of the backward endpoint,
      !then the (coord1, coord2, coord3) of the forward endpoint (backward/forward relative to direction of field)
      !The default is the spherical check, which compares r at both endpoints

          INTEGER :: idx_t
          REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
          BYTE :: connection_type
          REAL(num) :: Q_sign

          if ((pos_endpoints(1,idx_t) .lt. closed_fl_size)) then                 !one end at photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              connection_type = 0   !Closed
              Q_sign=1.0_num
            else                                                                 !second end away from photosphere
              connection_type = 1   !Open
              Q_sign=-1.0_num
            end if
          else                                                                   !one end away from photosphere
            if ((pos_endpoints(4,idx_t).lt. closed_fl_size)) then                !second end at photosphere
              connection_type = 1   !Open
              Q_sign=-1.0_num
            else                                                                 !second end away from photosphere
              connection_type = 2   !disconnected
              Q_sign=-1.0_num
            end if
          end if

      end subroutine fieldline_classification_user


      subroutine trace_cartesian_user(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                      B_gradB_interp,single_step,final_step,fl_classification, &
                                      pos_start,idx_t,pos_endpoints,pos_Q,pos_fieldline, &
                                      pos_step_start,pos_step_total,dl)
      !Trace fieldlines, calculate user defined quantity (Cartesian)

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        procedure(single_step_iface) :: final_step
        procedure(fl_classification_iface) :: fl_classification
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
          call final_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
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
          call final_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        call Bfull_interp(idx_in, pos_next, B_out)
        B_squared_end2 = B_out(1)**2+B_out(2)**2+B_out(3)**2

        fieldline_user(1,idx_t) = B_squared_end1/B_squared_end2

      end subroutine trace_cartesian_user


      subroutine trace_spherical_user(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                      B_gradB_interp,single_step,final_step,fl_classification, &
                                      pos_start,idx_t,pos_endpoints,pos_Q,pos_fieldline, &
                                      pos_step_start,pos_step_total,dl)
      !Trace fieldlines, calculate user defined quantity (spherical)

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        procedure(single_step_iface) :: final_step
        procedure(fl_classification_iface) :: fl_classification
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
          call final_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
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
          call final_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        call Bfull_interp(idx_in, pos_next, B_out)
        B_squared_end2 = B_out(1)**2+B_out(2)**2+B_out(3)**2

        fieldline_user(1,idx_t) = B_squared_end1/B_squared_end2

      end subroutine trace_spherical_user


      subroutine trace_cylindrical_user(check_position,intercept_boundary,B_interp,Bfull_interp, &
                                      B_gradB_interp,single_step,final_step,fl_classification, &
                                      pos_start,idx_t,pos_endpoints,pos_Q,pos_fieldline, &
                                      pos_step_start,pos_step_total,dl)
      !Trace fieldlines, calculate user defined quantity (cylindrical)

        procedure(check_position_iface) :: check_position
        procedure(intercept_boundary_iface) :: intercept_boundary
        procedure(B_interp_iface) :: B_interp
        procedure(Bfull_interp_iface) :: Bfull_interp
        procedure(B_gradB_interp_iface) :: B_gradB_interp
        procedure(single_step_iface) :: single_step
        procedure(single_step_iface) :: final_step
        procedure(fl_classification_iface) :: fl_classification
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
          call final_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
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
          call final_step(idx_in, pos_curr(:), pos_next(:), u_0, v_0, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        call Bfull_interp(idx_in, pos_next, B_out)
        B_squared_end2 = B_out(1)**2+B_out(2)**2+B_out(3)**2

        fieldline_user(1,idx_t) = B_squared_end1/B_squared_end2

      end subroutine trace_cylindrical_user


      subroutine cleanup_user_defined
      !Will be executed at the end; use this to deallocate, or do any other housekeeping on exit

        IF (ALLOCATED(fieldline_user)) DEALLOCATE(fieldline_user)

      end subroutine cleanup_user_defined


end module UFiT_User_Functions

