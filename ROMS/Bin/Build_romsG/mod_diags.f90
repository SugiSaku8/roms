      MODULE mod_diags
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  avgzeta   Time-averaged free surface                                !
!                                                                      !
!  Diagnostics fields for output.                                      !
!                                                                      !
!  DiaBio2d  Diagnostics for 2D biological terms.                      !
!  DiaBio3d  Diagnostics for 3D biological terms.                      !
!  DiaBio4d  Diagnostics for 4D bio-optical terms.                     !
!  DiaTrc    Diagnostics for tracer terms.                             !
!  DiaU2d    Diagnostics for 2D U-momentum terms.                      !
!  DiaV2d    Diagnostics for 2D V-momentum terms.                      !
!  DiaU3d    Diagnostics for 3D U-momentum terms.                      !
!  DiaV3d    Diagnostics for 3D V-momentum terms.                      !
!                                                                      !
!  Diagnostics fields work arrays.                                     !
!                                                                      !
!  DiaTwrk    Diagnostics work array for tracer terms.                 !
!  DiaU2wrk   Diagnostics work array for 2D U-momentum terms.          !
!  DiaV2wrk   Diagnostics work array for 2D V-momentum terms.          !
!  DiaRUbar   Diagnostics RHS array for 2D U-momentum terms.           !
!  DiaRVbar   Diagnostics RHS array for 2D V-momentum terms.           !
!  DiaU2int   Diagnostics array for 2D U-momentum terms                !
!               integrated over short barotropic timesteps.            !
!  DiaV2int   Diagnostics array for 2D U-momentum terms                !
!               integrated over short barotropic timesteps.            !
!  DiaRUfrc   Diagnostics forcing array for 2D U-momentum terms.       !
!  DiaRVfrc   Diagnostics forcing array for 2D V-momentum terms.       !
!  DiaU3wrk   Diagnostics work array for 3D U-momentum terms.          !
!  DiaV3wrk   Diagnostics work array for 3D V-momentum terms.          !
!  DiaRU      Diagnostics RHS array for 3D U-momentum terms.           !
!  DiaRV      Diagnostics RHS array for 3D V-momentum terms.           !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
!
        implicit none
!
        PUBLIC :: allocate_diags
        PUBLIC :: deallocate_diags
        PUBLIC :: initialize_diags
!
!-----------------------------------------------------------------------
!  Define T_DIAGS structure.
!-----------------------------------------------------------------------
!
        TYPE T_DIAGS
          real(r8), pointer :: avgzeta(:,:)
        END TYPE T_DIAGS
!
        TYPE (T_DIAGS), allocatable :: DIAGS(:)
!
      CONTAINS
!
      SUBROUTINE allocate_diags (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1 ) allocate ( DIAGS(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Diagnostic arrays.
!
      allocate ( DIAGS(ng) % avgzeta(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
!
      RETURN
      END SUBROUTINE allocate_diags
!
      SUBROUTINE deallocate_diags (ng)
!
!=======================================================================
!                                                                      !
!  This routine deallocates all variables in the module for all nested !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Modules/mod_diags.F"//", deallocate_diags"
!
!-----------------------------------------------------------------------
!  Deallocate derived-type DIAGS structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(DIAGS)) deallocate ( DIAGS )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_diags
!
      SUBROUTINE initialize_diags (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, idiag, j
      integer :: itrc, iband, k
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!  Set array initialization range.
!
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          DIAGS(ng) % avgzeta(i,j) = IniVal
        END DO
      END DO
!
      RETURN
      END SUBROUTINE initialize_diags
      END MODULE mod_diags
