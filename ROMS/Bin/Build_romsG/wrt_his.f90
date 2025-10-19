      MODULE wrt_his_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module writes requested model fields into the HISTORY output   !
!  file using either the standard NetCDF library or the Parallel-IO    !
!  (PIO) library.                                                      !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE extract_slice_mod,     ONLY : extract_slice
      USE nf_fwrite2d_mod,       ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod,       ONLY : nf_fwrite3d
      USE omega_mod,             ONLY : scale_omega
      USE strings_mod,           ONLY : FoundError
      USE uv_rotate_mod,         ONLY : uv_rotate2d
      USE uv_rotate_mod,         ONLY : uv_rotate3d
      USE vorticity_mod,         ONLY : pvorticity3d, rvorticity3d
!
      implicit none
!
      PUBLIC  :: wrt_his
      PRIVATE :: wrt_his_nf90
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_his (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: IminS, ImaxS, JminS, JmaxS
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_his.F"
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
      SELECT CASE (HIS(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_his_nf90 (ng, iNLM, tile,                            &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS)
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) HIS(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, 167, MyFile)) RETURN
!
  10  FORMAT (' WRT_HIS - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_his
!
!***********************************************************************
      SUBROUTINE wrt_his_nf90 (ng, model, tile,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Fcount, gfactor, gtype, ifield, status
      integer :: i, itrc, j, k
!
      real(dp) :: scale
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
      real(r8), allocatable :: Fr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_his.F"//", wrt_his_nf90"
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
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 222, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set time record index.
!
      HIS(ng)%Rindex=HIS(ng)%Rindex+1
      Fcount=HIS(ng)%load
      HIS(ng)%Nrec(Fcount)=HIS(ng)%Nrec(Fcount)+1
!
!  Report.
!
      IF (Master) WRITE (stdout,10) kstp(ng), nrhs(ng), HIS(ng)%Rindex
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, model, HIS(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/HIS(ng)%Rindex/), (/1/),                  &
     &                      ncid = HIS(ng)%ncid,                        &
     &                      varid = HIS(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, 262, MyFile)) RETURN
!
!  Write time-varying depths of RHO-points.
!
      IF (Hout(idpthR,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthR,             &
     &                     HIS(ng)%Vid(idpthR),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % z_r,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 370, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthR)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of U-points.
!
      IF (Hout(idpthU,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthU,             &
     &                     HIS(ng)%Vid(idpthU),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 402, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthU)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of V-points.
!
      IF (Hout(idpthV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthV,             &
     &                     HIS(ng)%Vid(idpthV),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask,                            &
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 434, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthV)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of W-points.
!
      IF (Hout(idpthW,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthW,             &
     &                     HIS(ng)%Vid(idpthW),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % z_w,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 458, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthW)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idFsur,             &
     &                     HIS(ng)%Vid(idFsur),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % zeta(:,:,kstp(ng)))
        IF (FoundError(status, nf90_noerr, 487, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUbar,             &
     &                     HIS(ng)%Vid(idUbar),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_full,                       &
     &                     OCEAN(ng) % ubar(:,:,kstp(ng)))
        IF (FoundError(status, nf90_noerr, 552, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVbar,             &
     &                     HIS(ng)%Vid(idVbar),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_full,                       &
     &                     OCEAN(ng) % vbar(:,:,kstp(ng)))
        IF (FoundError(status, nf90_noerr, 671, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_full,                        &
     &                    OCEAN(ng) % ubar(:,:,kstp(ng)),               &
     &                    OCEAN(ng) % vbar(:,:,kstp(ng)),               &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idu2dE,             &
     &                     HIS(ng)%Vid(idu2dE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_full,                       &
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, 810, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idv2dN,             &
     &                     HIS(ng)%Vid(idv2dN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_full,                       &
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, 827, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idUvel,             &
     &                     HIS(ng)%Vid(idUvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask_full,                       &
     &                     OCEAN(ng) % u(:,:,:,nrhs(ng)))
        IF (FoundError(status, nf90_noerr, 854, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D U-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idUzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*u3dvar
!
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % u(:,:,:,nrhs(ng)),              &
     &                      GRID(ng) % z_v,                             &
     &                      GRID(ng) % umask_full,                      &
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idUzsl,             &
     &                     HIS(ng)%Vid(idUzsl),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % umask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 918, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVvel,             &
     &                     HIS(ng)%Vid(idVvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask_full,                       &
     &                     OCEAN(ng) % v(:,:,:,nrhs(ng)))
        IF (FoundError(status, nf90_noerr, 967, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idVzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*v3dvar
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % v(:,:,:,nrhs(ng)),              &
     &                      GRID(ng) % z_v,                             &
     &                      GRID(ng) % vmask_full,                      &
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVzsl,             &
     &                     HIS(ng)%Vid(idVzsl),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % vmask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 1031, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid.
!
      IF (Hout(idu3dE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idu3dE,             &
     &                     HIS(ng)%Vid(idu3dE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_full,                       &
     &                     OCEAN(ng) % ua)
        IF (FoundError(status, nf90_noerr, 1080, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idUzsE,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*r3dvar
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % ua,                             &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % rmask,                           &
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idUzsE,             &
     &                     HIS(ng)%Vid(idUzsE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % rmask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 1117, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid.
!
      IF (Hout(idv3dN,ng)) THEN
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idv3dN,             &
     &                     HIS(ng)%Vid(idv3dN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_full,                       &
     &                     OCEAN(ng) % va)
        IF (FoundError(status, nf90_noerr, 1139, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idVzsN,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*r3dvar
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % va,                             &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % rmask,                           &
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVzsN,             &
     &                     HIS(ng)%Vid(idVzsN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % rmask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 1176, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out potential vorticity (m-1 s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dPV,ng).or.(Hout(idPVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*p3dvar
!
        CALL pvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS, nrhs(ng),        &
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % umask,   GRID(ng) % vmask,        &
                           GRID(ng) % f,                                &
                           GRID(ng) % om_u,    GRID(ng) % on_v,         &
                           GRID(ng) % pm,      GRID(ng) % pn,           &
                           GRID(ng) % z_r,     OCEAN(ng) % pden,        &
                           OCEAN(ng) % u,      OCEAN(ng) % v,           &
                           Fr3d)
!
        IF (Hout(id3dPV,ng)) THEN
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, id3dPV,           &
     &                       HIS(ng)%Vid(id3dPV),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % pmask,                          &
     &                       Fr3d)
          IF (FoundError(status, nf90_noerr, 1221, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dPV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idPVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                  &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fr3d,                                     &
     &                        GRID(ng) % z_v,                           &
     &                        GRID(ng) % pmask,                         &
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idPVzs,           &
     &                       HIS(ng)%Vid(idPVzs),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
     &                       GRID(ng) % pmask,                          &
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, 1265, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idPVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out relative vorticity (s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dRV,ng).or.(Hout(idRVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*p3dvar
!
        CALL rvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj, nrhs(ng),                &
     &                     GRID(ng) % pmask,                            &
                           GRID(ng) % om_u, GRID(ng) % on_v,            &
                           GRID(ng) % pm,   GRID(ng) % pn,              &
                           OCEAN(ng) % u,   OCEAN(ng) % v,              &
                           Fr3d)
!
        IF (Hout(id3dRV,ng)) THEN
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, id3dRV,           &
     &                       HIS(ng)%Vid(id3dRV),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % pmask,                          &
     &                       Fr3d)
          IF (FoundError(status, nf90_noerr, 1308, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dRV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idRVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                  &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fr3d,                                     &
     &                        GRID(ng) % z_v,                           &
     &                        GRID(ng) % pmask,                         &
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idRVzs,           &
     &                       HIS(ng)%Vid(idRVzs),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
     &                       GRID(ng) % pmask,                          &
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, 1352, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idRVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*w3dvar
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d)
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idOvel,             &
     &                     HIS(ng)%Vid(idOvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 1387, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idWvel,             &
     &                     HIS(ng)%Vid(idWvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % wvel)
        IF (FoundError(status, nf90_noerr, 1447, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idWvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idTvar(itrc),     &
     &                       HIS(ng)%Tid(itrc),                         &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % rmask,                          &
     &                       OCEAN(ng) % t(:,:,:,nrhs(ng),itrc))
          IF (FoundError(status, nf90_noerr, 1471, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out tracer type variables at specified constant depth slices.
!
      DO itrc=1,NT(ng)
        IF (Hout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
            Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
          END IF
          scale=1.0_dp
          gtype=gfactor*r3dvar
          CALL extract_slice (ng, model, tile, gtype,                   &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng) % t(:,:,:,nrhs(ng),itrc),       &
     &                        GRID(ng) % z_r,                           &
     &                        GRID(ng) % rmask,                         &
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idzslT(itrc),     &
     &                       HIS(ng)%Vid(idzslT(itrc)),                 &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
     &                       GRID(ng) % rmask,                          &
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, 1510, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idzslT(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          IF (itrc.eq.NT(ng)) deallocate (Wr3d)
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idDano,             &
     &                     HIS(ng)%Vid(idDano),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % rho)
        IF (FoundError(status, nf90_noerr, 1564, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVvis,             &
     &                     HIS(ng)%Vid(idVvis),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1665, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idTdif,             &
     &                     HIS(ng)%Vid(idTdif),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1689, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface active tracers fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp                   ! Celsius m/s to W/m2
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_dp
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idTsur(itrc),     &
     &                       HIS(ng)%Vid(idTsur(itrc)),                 &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (FoundError(status, nf90_noerr, 2016, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idEmPf,             &
     &                     HIS(ng)%Vid(idEmPf),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % stflux(:,:,isalt))
        IF (FoundError(status, nf90_noerr, 2164, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idEmPf)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=rho0                          ! m2/s2 to Pa
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUsms,             &
     &                     HIS(ng)%Vid(idUsms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     FORCES(ng) % sustr)
        IF (FoundError(status, nf90_noerr, 2217, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVsms,             &
     &                     HIS(ng)%Vid(idVsms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     FORCES(ng) % svstr)
        IF (FoundError(status, nf90_noerr, 2244, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUbms,             &
     &                     HIS(ng)%Vid(idUbms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     FORCES(ng) % bustr)
        IF (FoundError(status, nf90_noerr, 2267, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVbms,             &
     &                     HIS(ng)%Vid(idVbms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     FORCES(ng) % bvstr)
        IF (FoundError(status, nf90_noerr, 2290, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, HIS(ng)%name, HIS(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 2366, MyFile)) RETURN
!
  10  FORMAT (2x,'WRT_HIS_NF90     - writing history', t42,             &
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
  20  FORMAT (/,' WRT_HIS_NF90 - error while writing variable: ',a,     &
     &        /,16x,'into history NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_his_nf90
!
      END MODULE wrt_his_mod
