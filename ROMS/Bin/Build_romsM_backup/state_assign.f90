      MODULE state_assign_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module associates/disassociates temporary work state variables !
!  "my_var" from ROMS kernel variables for efficient use of computer   !
!  memory.                                                             !
!                                                                      !
!  In the 4D-Var algorithms, several state/control vector copies are   !
!  needed to store ROMS trajectories, background, error covariance     !
!  standard deviations, error covariance normalization coefficients,   !
!  impulse forcing, eigenvectors, conjugate directions, and other      !
!  transformations. They are also used to read/write state vectors     !
!  into NetCDF files. The list below documents how these state vectors !
!  are used to help manage memory better in some parts of code.        !
!                                                                      !
!  STATE     ID    INDEX                                               !
!  --------  ----  -----                                               !
!  nl_state  iNLM  1,2     basic nonlinear kernel/prior state          !
!  ad_state  iADM  1,2     basic adjoint kernel/control state          !
!  tl_state  iTLM  1,2     basic tangent linear kernel/control state   !
!            iRPM  1,2     basic representer state                     !
!  b_state   14    1       initial conditions B normalization          !
!            15    2       model error B normalization                 !
!  e_state   10    1       initial conditions B standard deviations    !
!            11    1       model error B standard deviation            !
!  d_state   18    N/A     conjugate direction, eigenvectors           !
!  f_state   7     N/A     time interpolated impulse forcing           !
!  fG_state  19    1,2     impulse forcing snapshots                   !
!  fS_state  N/A   1,2     time convolutions snapshots (not used here) !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ocean
!
      implicit none
!
      PUBLIC :: state_assign
      PUBLIC :: state_unassign
!
!  Declare pointers for temporary work state variables.
!
      real (r8), pointer :: my_zeta(:,:)  => NULL()
      real (r8), pointer :: my_ubar(:,:)  => NULL()
      real (r8), pointer :: my_vbar(:,:)  => NULL()
      real (r8), pointer :: my_t(:,:,:,:) => NULL()
      real (r8), pointer :: my_u(:,:,:)   => NULL()
      real (r8), pointer :: my_v(:,:,:)   => NULL()
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CONTAINS
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     
!***********************************************************************
      SUBROUTINE state_assign (ng, Istate, TimeIndex)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Istate, TimeIndex
!
!-----------------------------------------------------------------------
!  Assign temporarily work state variables to requested ROMS state
!  vector.
!-----------------------------------------------------------------------
!
!  Associate to required state vector variables.
!
      SELECT CASE (Istate)
        CASE (iNLM)
          IF (associated(OCEAN(ng)%zeta)) THEN
            my_zeta => OCEAN(ng) % zeta(:,:,TimeIndex)
            my_ubar => OCEAN(ng) % ubar(:,:,TimeIndex)
            my_vbar => OCEAN(ng) % vbar(:,:,TimeIndex)
            my_t    => OCEAN(ng) % t(:,:,:,TimeIndex,:)
            my_u    => OCEAN(ng) % u(:,:,:,TimeIndex)
            my_v    => OCEAN(ng) % v(:,:,:,TimeIndex)
            IF (Master) WRITE (stdout,10) 'Nonlinear', TimeIndex
          ELSE
            IF (Master) WRITE (stdout,30) 'Nonlinear'
          END IF
      END SELECT
!
 10   FORMAT (/,2x,'STATE_ASSIGN     - Work state variables ',          &
     &        'pointers are associated with => ',a,' variables, ',      &
     &        'TimeIndex = ',i0)
 20   FORMAT (/,2x,'STATE_ASSIGN     - Work state variables ',          &
     &        'pointers are associated with => ',a,' variables.')
 30   FORMAT (/,2x,'STATE_ASSIGN     - unallocated ',a,                 &
     &        'state variables.')
!
      RETURN
      END SUBROUTINE state_assign
!     
!***********************************************************************
      SUBROUTINE state_unassign
!***********************************************************************
!
!-----------------------------------------------------------------------
!  Nullify temporary work state variables.
!-----------------------------------------------------------------------
!
      IF (associated(my_zeta)) nullify (my_zeta)
      IF (associated(my_ubar)) nullify (my_ubar)
      IF (associated(my_vbar)) nullify (my_vbar)
      IF (associated(my_t))    nullify (my_t)
      IF (associated(my_u))    nullify (my_u)
      IF (associated(my_v))    nullify (my_v)
      IF (Master) WRITE(stdout,10)
 10   FORMAT (/,2x,'STATE_UNASSIGN   - Work state variables ',          &
     &        'pointers nullified.',/)
!
      RETURN
      END SUBROUTINE state_unassign
!
      END MODULE state_assign_mod
