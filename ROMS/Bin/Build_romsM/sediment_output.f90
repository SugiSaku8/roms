      MODULE sediment_output_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module defines/writes Waves Effect on Currents variables into  !
!  output NetCDF files.                                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_average
      USE mod_bbl
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_sedbed
      USE mod_sediment
      USE mod_stepping
!
      USE def_var_mod,     ONLY : def_var
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
      USE omega_mod,       ONLY : scale_omega
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC :: sediment_def_nf90
      PUBLIC :: sediment_wrt_nf90
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE sediment_def_nf90 (ng, model, ldef, VarOut, S,         &
     &                              t2dgrd, u2dgrd, v2dgrd,             &
     &                              b3dgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: t2dgrd(:), u2dgrd(:), v2dgrd(:)
      integer, intent(in), optional :: b3dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25
      integer :: i, itrc, j, nvd3, nvd4, status
!
      real(r8) :: Aval(6)
!
      character (len=13)  :: Prefix
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/Sediment/sediment_output.F"//", sediment_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sediment output variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 123, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
!
!  Set long name prefix string.
!
        Prefix=CHAR(32)                                   ! blank space
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!  Define exposed sediment layer properties.
!
        DO i=1,MBOTP
          IF (VarOut(idBott(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBott(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idBott(i)))
            ELSE
              Vinfo( 2)=Vname(2,idBott(i))
            END IF
            Vinfo( 3)=Vname(3,idBott(i))
            Vinfo(14)=Vname(4,idBott(i))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idBott(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idBott(i),ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idBott(i)), NF_FOUT,               &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 348, MyFile)) RETURN
          END IF
        END DO
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize local logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Waves Effect on Currents variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
          END IF
          DO itrc=1,MBOTP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              S(ng)%Vid(idBott(itrc))=var_id(i)
            END IF
          END DO
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO i=1,MBOTP
          IF (.not.got_var(idBott(i)).and.VarOut(idBott(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idBott(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
      END IF QUERY
!
  10  FORMAT (1pe11.4,1x,'millimeter')
  20  FORMAT (/,' SEDIMENT_DEF_NF90 - unable to find variable: ',       &
     &        a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE sediment_def_nf90
!
!***********************************************************************
      SUBROUTINE sediment_wrt_nf90 (ng, model, tile,                    &
    &                               LBi, UBi, LBj, UBj,                 &
    &                               VarOut, S)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Linstataneous
!
      integer :: gfactor, gtype, i, status
!
      real(dp) :: scale
!
      character (len=256) :: ncname
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/Sediment/sediment_output.F"//", sediment_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Currents output variables into specified
!  output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 1031, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set instantaneous fields.
!
      IF ((S(ng)%ncid.eq.S(ng)%ncid).or.                                &
     &    (S(ng)%ncid.eq.QCK(ng)%ncid)) THEN
        Linstataneous=.TRUE.
      ELSE
        Linstataneous=.FALSE.                  ! time-averged fiels
      END IF
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBOTP
        IF (VarOut(idBott(i),ng)) THEN
          IF (i.eq.itauc) THEN
            scale=rho0
          ELSE
            scale=1.0_dp
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idBott(i),          &
     &                       S(ng)%Vid(idBott(i)),                      &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       SEDBED(ng) % bottom(:,:,i))
          IF (FoundError(status, nf90_noerr, 1262, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBott(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
  10  FORMAT (/," SEDIMENT_WRT_NF90 - error while writing variable '",  &
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE sediment_wrt_nf90
!
      END MODULE sediment_output_mod
