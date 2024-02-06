

module UFiT_User_Functions
      USE UFiT_Definitions_Fortran
      implicit none


      contains


      subroutine prepare_user_defined
      !Will be executed if running in user mode

        ALLOCATE(fieldline_user(1,numin_tot))

      end subroutine prepare_user_defined


      subroutine trace_cartesian_user(check_position,intercept_boundary,B_interp,B_gradB_interp, &
                                 single_step,pos_start,idx_t,pos_endpoints,pos_Q,pos_fieldline, &
                                 pos_step_start,pos_step_total,dl)
      !Trace fieldlines, calculate user defined quantity (Cartesian)

        interface
          subroutine check_position (pos_in,pos_out,res)
            IMPORT :: num
            REAL(num) :: pos_in(3), pos_out(3)
            LOGICAL :: res
          end subroutine check_position
        end interface
        interface
          subroutine intercept_boundary (pos_in,direction,delta_out)
            IMPORT :: num
            REAL(num) :: pos_in(3), direction(3), delta_out
          end subroutine intercept_boundary
        end interface
        interface
          subroutine B_interp (idx1, idx2, idx3, pos_in, B_out)
            IMPORT :: num
            INTEGER :: idx1, idx2, idx3
            REAL(num), INTENT(IN) :: pos_in(3)
            REAL(num), INTENT(OUT) :: B_out(3)
          end subroutine B_interp
        end interface
        interface
          subroutine B_gradB_interp (idx1, idx2, idx3, pos_in, B_out, gradB_out)
            IMPORT :: num
            INTEGER :: idx1, idx2, idx3
            REAL(num), INTENT(IN) :: pos_in(3)
            REAL(num), INTENT(OUT) :: B_out(3)
            REAL(num), INTENT(OUT) :: gradB_out(3,3)
          end subroutine B_gradB_interp
        end interface
        interface
          subroutine single_step (idx1, idx2, idx3, pos_in, pos_out, u_vec, v_vec, dl, &
                                  keep_running, check_position, B_interp, B_gradB_interp)
            IMPORT :: num
            INTEGER :: idx1, idx2, idx3
            REAL(num) :: pos_in(3), pos_out(3), u_vec(3), v_vec(3), dl
            LOGICAL :: keep_running
            interface
              subroutine check_position (pos_in,pos_out,res)
                IMPORT :: num
                REAL(num) :: pos_in(3), pos_out(3)
                LOGICAL :: res
              end subroutine check_position
            end interface
            interface
              subroutine B_interp (idx1, idx2, idx3, pos_in, B_out)
                IMPORT :: num
                INTEGER :: idx1, idx2, idx3
                REAL(num), INTENT(IN) :: pos_in(3)
                REAL(num), INTENT(OUT) :: B_out(3)
              end subroutine B_interp
            end interface
            interface
              subroutine B_gradB_interp (idx1, idx2, idx3, pos_in, B_out, gradB_out)
                IMPORT :: num
                INTEGER :: idx1, idx2, idx3
                REAL(num), INTENT(IN) :: pos_in(3)
                REAL(num), INTENT(OUT) :: B_out(3)
                REAL(num), INTENT(OUT) :: gradB_out(3,3)
              end subroutine B_gradB_interp
            end interface
          end subroutine single_step
        end interface

        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_x, idx_y, idx_z
        LOGICAL :: keep_running
        REAL(num) :: pos_out(3), B_out(3), mod_Bout, last_stepsize
        REAL(num) :: pos_next(3), pos_curr(3), u_0(3), v_0(3)
        REAL(num) :: B_squared_end1, B_squared_end2

        idx_x = 1
        idx_y = 1
        idx_z = 1

      !Backwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_x, idx_y, idx_z, pos_curr(:), pos_next(:), u_0, v_0, -dl, &
                               keep_running, check_position, B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_x, idx_y, idx_z, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_x, idx_y, idx_z, pos_curr(:), pos_next(:), u_0, v_0, &
                              -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                              B_interp, B_gradB_interp)
        end if
        call B_interp(idx_x, idx_y, idx_z, pos_next, B_out)
        B_squared_end1 = B_out(1)**2+B_out(2)**2+B_out(3)**2

      !Forwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_x, idx_y, idx_z, pos_curr(:), pos_next(:), u_0, v_0, &
                               dl, keep_running, check_position, B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_x, idx_y, idx_z, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_x, idx_y,idx_z,pos_curr(:), pos_next(:), u_0, v_0, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        call B_interp(idx_x, idx_y, idx_z, pos_next, B_out)
        B_squared_end2 = B_out(1)**2+B_out(2)**2+B_out(3)**2

        fieldline_user(1,numin_tot) = B_squared_end1/B_squared_end2

      end subroutine trace_cartesian_user


      subroutine trace_spherical_user(check_position,intercept_boundary,B_interp,B_gradB_interp, &
                                      single_step,pos_start,idx_t,pos_endpoints,pos_Q, &
                                      pos_fieldline,pos_step_start,pos_step_total,dl)
      !Trace fieldlines, calculate user defined quantity (spherical)

        interface
          subroutine check_position (pos_in,pos_out,res)
            IMPORT :: num
            REAL(num) :: pos_in(3), pos_out(3)
            LOGICAL :: res
          end subroutine check_position
        end interface
        interface
          subroutine intercept_boundary (pos_in,direction,delta_out)
            IMPORT :: num
            REAL(num) :: pos_in(3), direction(3), delta_out
          end subroutine intercept_boundary
        end interface
        interface
          subroutine B_interp (idx1, idx2, idx3, pos_in, B_out)
            IMPORT :: num
            INTEGER :: idx1, idx2, idx3
            REAL(num), INTENT(IN) :: pos_in(3)
            REAL(num), INTENT(OUT) :: B_out(3)
          end subroutine B_interp
        end interface
        interface
          subroutine B_gradB_interp (idx1, idx2, idx3, pos_in, B_out, gradB_out)
            IMPORT :: num
            INTEGER :: idx1, idx2, idx3
            REAL(num), INTENT(IN) :: pos_in(3)
            REAL(num), INTENT(OUT) :: B_out(3)
            REAL(num), INTENT(OUT) :: gradB_out(3,3)
          end subroutine B_gradB_interp
        end interface
        interface
          subroutine single_step (idx1, idx2, idx3, pos_in, pos_out, u_vec, v_vec, dl, &
                                  keep_running, check_position, B_interp, B_gradB_interp)
            IMPORT :: num
            INTEGER :: idx1, idx2, idx3
            REAL(num) :: pos_in(3), pos_out(3), u_vec(3), v_vec(3), dl
            LOGICAL :: keep_running
            interface
              subroutine check_position (pos_in,pos_out,res)
                IMPORT :: num
                REAL(num) :: pos_in(3), pos_out(3)
                LOGICAL :: res
              end subroutine check_position
            end interface
            interface
              subroutine B_interp (idx1, idx2, idx3, pos_in, B_out)
                IMPORT :: num
                INTEGER :: idx1, idx2, idx3
                REAL(num), INTENT(IN) :: pos_in(3)
                REAL(num), INTENT(OUT) :: B_out(3)
              end subroutine B_interp
            end interface
            interface
              subroutine B_gradB_interp (idx1, idx2, idx3, pos_in, B_out, gradB_out)
                IMPORT :: num
                INTEGER :: idx1, idx2, idx3
                REAL(num), INTENT(IN) :: pos_in(3)
                REAL(num), INTENT(OUT) :: B_out(3)
                REAL(num), INTENT(OUT) :: gradB_out(3,3)
              end subroutine B_gradB_interp
            end interface
          end subroutine single_step
        end interface

        REAL(num) :: pos_start(3)
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: pos_endpoints
        REAL(num), DIMENSION(:), ALLOCATABLE :: pos_Q
        REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: pos_fieldline
        INTEGER :: idx_t
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_start
        INTEGER, DIMENSION(:), ALLOCATABLE :: pos_step_total
        REAL(num) :: dl

        INTEGER :: step_num, idx_r, idx_th, idx_p
        LOGICAL :: keep_running
        REAL(num) :: pos_out(3), B_out(3), mod_Bout, last_stepsize
        REAL(num) :: pos_next(3), pos_curr(3), u_0(3), v_0(3)
        REAL(num) :: B_squared_end1, B_squared_end2

        idx_r = 1
        idx_th = 1
        idx_p = 1

      !Backwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_r, idx_th, idx_p, pos_curr(:), pos_next(:), u_0, v_0, -dl, &
                               keep_running, check_position,B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_r, idx_th, idx_p, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = -B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_r, idx_th, idx_p, pos_curr(:), pos_next(:), u_0, v_0, &
                              -last_stepsize*dl*0.9999_num, keep_running, check_position, &
                              B_interp, B_gradB_interp)
        end if
        call B_interp(idx_r, idx_th, idx_p, pos_next, B_out)
        B_squared_end1 = B_out(1)**2+B_out(2)**2+B_out(3)**2

      !Forwards integration
        pos_next(:) = pos_start(:)
        step_num = 1
        keep_running = .true.
        do while (keep_running .and. (step_num .le. MAX_STEPS))
          pos_curr(:) = pos_next(:)
          call single_step(idx_r, idx_th, idx_p, pos_curr(:), pos_next(:), u_0, v_0, &
                               dl, keep_running, check_position,B_interp, B_gradB_interp)
          step_num = step_num + 1
        end do
        if (.not. keep_running) then
          call check_position(pos_curr(:),pos_out,keep_running)
          call B_interp(idx_r, idx_th, idx_p, pos_out, B_out)
          mod_Bout = SQRT(B_out(1)**2+B_out(2)**2+B_out(3)**2)
          B_out(:) = B_out(:)*dl/mod_Bout
          call intercept_boundary(pos_curr(:),B_out,last_stepsize)
          call single_step(idx_r, idx_th, idx_p, pos_curr(:), pos_next(:), u_0, v_0, &
                               last_stepsize*dl*0.9999_num, keep_running, check_position, &
                               B_interp, B_gradB_interp)
        end if
        call B_interp(idx_r, idx_th, idx_p, pos_next, B_out)
        B_squared_end2 = B_out(1)**2+B_out(2)**2+B_out(3)**2

        fieldline_user(1,numin_tot) = B_squared_end1/B_squared_end2

      end subroutine trace_spherical_user


end module UFiT_User_Functions

