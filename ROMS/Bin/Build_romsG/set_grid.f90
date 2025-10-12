      SUBROUTINE set_grid (ng, model)
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine sets application grid and associated variables and     !
!  parameters. It called only once during the initialization stage.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
      USE analytical_mod
      USE distribute_mod,       ONLY : mp_bcasti
      USE get_grid_mod,         ONLY : get_grid
      USE get_nudgcoef_mod,     ONLY : get_nudgcoef
      USE metrics_mod,          ONLY : metrics
      USE strings_mod,          ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: tile
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/set_grid.F"
!
!-----------------------------------------------------------------------
!  Set horizontal grid, bathymetry, and Land/Sea masking (if any).
!  Use analytical functions or read in from a grid NetCDF.
!-----------------------------------------------------------------------
!
      CALL get_grid (ng, MyRank, model)
      CALL mp_bcasti (ng, model, exit_flag)
      IF (FoundError(exit_flag, NoError, 88, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Set vertical terrain-following coordinate transformation function.
!-----------------------------------------------------------------------
!
      CALL set_scoord (ng)
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
      CALL set_weights (ng)
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL metrics (ng, tile, model)
      END DO
!
!-----------------------------------------------------------------------
!  If appropriate, set spatially varying nudging coefficients time
!  scales.
!-----------------------------------------------------------------------
!
      IF (Lnudging(ng)) THEN
        CALL get_nudgcoef (ng, MyRank, model)
        CALL mp_bcasti (ng, model, exit_flag)
        IF (FoundError(exit_flag, NoError, 190, MyFile)) RETURN
      END IF
!
      RETURN
      END SUBROUTINE set_grid
