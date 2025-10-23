      MODULE bbl_output_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module defines/writes Bottom Boundary Layer model variables    !
!  into output NetCDF files.                                           !
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
      PUBLIC :: bbl_def_nf90
      PUBLIC :: bbl_wrt_nf90
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE bbl_def_nf90 (ng, model, ldef, VarOut, S,              &
     &                         t2dgrd, u2dgrd, v2dgrd)
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
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25
      integer :: i, j, nvd3, nvd4, status
!
      real(r8) :: Aval(6)
!
      character (len=13)  :: Prefix
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/BBL/bbl_output.F"//", bbl_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define Bottom Boundary Layer (BBL) and Waves output variables.
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
!  Define wind-induced bottom orbital velocity.
!
        IF (VarOut(idWorb,ng)) THEN
          Vinfo( 1)=Vname(1,idWorb)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWorb))
          ELSE
            Vinfo( 2)=Vname(2,idWorb)
          END IF
          Vinfo( 3)=Vname(3,idWorb)
          Vinfo(14)=Vname(4,idWorb)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idWorb)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWorb,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWorb),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 182, MyFile)) RETURN
        END IF
!
!  Define bottom U-current stress.
!
        IF (VarOut(idUbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbrs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbrs))
          ELSE
            Vinfo( 2)=Vname(2,idUbrs)
          END IF
          Vinfo( 3)=Vname(3,idUbrs)
          Vinfo(14)=Vname(4,idUbrs)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbrs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbrs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbrs),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 208, MyFile)) RETURN
        END IF
!
!  Define bottom V-current stress.
!
        IF (VarOut(idVbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbrs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbrs))
          ELSE
            Vinfo( 2)=Vname(2,idVbrs)
          END IF
          Vinfo( 3)=Vname(3,idVbrs)
          Vinfo(14)=Vname(4,idVbrs)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbrs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbrs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbrs),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 231, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bottom U-wave stress.
!
        IF (VarOut(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbws))
          ELSE
            Vinfo( 2)=Vname(2,idUbws)
          END IF
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbws)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbws,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbws),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 254, MyFile)) RETURN
        END IF
!
!  Define bottom wind-induced, bottom V-wave stress.
!
        IF (VarOut(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbws))
          ELSE
            Vinfo( 2)=Vname(2,idVbws)
          END IF
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbws)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbws,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbws),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 277, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom U-wave stress.
!
        IF (VarOut(idUbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbcs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbcs))
          ELSE
            Vinfo( 2)=Vname(2,idUbcs)
          END IF
          Vinfo( 3)=Vname(3,idUbcs)
          Vinfo(14)=Vname(4,idUbcs)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbcs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbcs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbcs),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 300, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom V-wave stress.
!
        IF (VarOut(idVbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbcs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbcs))
          ELSE
            Vinfo( 2)=Vname(2,idVbcs)
          END IF
          Vinfo( 3)=Vname(3,idVbcs)
          Vinfo(14)=Vname(4,idVbcs)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbcs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbcs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbcs),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 323, MyFile)) RETURN
        END IF
!
!  Define maximum wave and current bottom stress magnitude.
!
        IF (VarOut(idUVwc,ng)) THEN
          Vinfo( 1)=Vname(1,idUVwc)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUVwc))
          ELSE
            Vinfo( 2)=Vname(2,idUVwc)
          END IF
          Vinfo( 3)=Vname(3,idUVwc)
          Vinfo(14)=Vname(4,idUVwc)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUVwc)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUVwc,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUVwc),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 346, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital U-velocity.
!
        IF (VarOut(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbot))
          ELSE
            Vinfo( 2)=Vname(2,idUbot)
          END IF
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbot)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbot,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbot),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 369, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital V-velocity.
!
        IF (VarOut(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbot))
          ELSE
            Vinfo( 2)=Vname(2,idVbot)
          END IF
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbot)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbot,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbot),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 392, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum above bed.
!
        IF (VarOut(idUbur,ng)) THEN
          Vinfo( 1)=Vname(1,idUbur)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbur))
          ELSE
            Vinfo( 2)=Vname(2,idUbur)
          END IF
          Vinfo( 3)=Vname(3,idUbur)
          Vinfo(14)=Vname(4,idUbur)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbur,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbur),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 415, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum above bed.
!
        IF (VarOut(idVbvr,ng)) THEN
          Vinfo( 1)=Vname(1,idVbvr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbvr))
          ELSE
            Vinfo( 2)=Vname(2,idVbvr)
          END IF
          Vinfo( 3)=Vname(3,idVbvr)
          Vinfo(14)=Vname(4,idVbvr)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbvr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbvr,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbvr),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 438, MyFile)) RETURN
        END IF
!
!  Define wind-induced mean wave direction.
!
        IF (VarOut(idWdir,ng)) THEN
          Vinfo( 1)=Vname(1,idWdir)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdir))
          ELSE
            Vinfo( 2)=Vname(2,idWdir)
          END IF
          Vinfo( 3)=Vname(3,idWdir)
          Vinfo(14)=Vname(4,idWdir)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idWdir)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdir,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdir),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 609, MyFile)) RETURN
        END IF
!
!  Define wind-induced peak wave direction.
!
        IF (VarOut(idWdip,ng)) THEN
          Vinfo( 1)=Vname(1,idWdip)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdip))
          ELSE
            Vinfo( 2)=Vname(2,idWdip)
          END IF
          Vinfo( 3)=Vname(3,idWdip)
          Vinfo(14)=Vname(4,idWdip)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idWdip)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdip,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdip),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 635, MyFile)) RETURN
        END IF
!
!  Define wind-induced bottom wave period.
!
        IF (VarOut(idWpbt,ng)) THEN
          Vinfo( 1)=Vname(1,idWpbt)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWpbt))
          ELSE
            Vinfo( 2)=Vname(2,idWpbt)
          END IF
          Vinfo( 3)=Vname(3,idWpbt)
          Vinfo(14)=Vname(4,idWpbt)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idWpbt)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWpbt,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWpbt),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 687, MyFile)) RETURN
        END IF
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
!  Bottom Boundary Layer model variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWorb))) THEN
            got_var(idWorb)=.TRUE.
            S(ng)%Vid(idWorb)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            S(ng)%Vid(idUbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            S(ng)%Vid(idVbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            S(ng)%Vid(idUbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            S(ng)%Vid(idVbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            S(ng)%Vid(idUbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            S(ng)%Vid(idVbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVwc))) THEN
            got_var(idUVwc)=.TRUE.
            S(ng)%Vid(idUVwc)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            S(ng)%Vid(idUbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            S(ng)%Vid(idVbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            S(ng)%Vid(idUbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            S(ng)%Vid(idVbvr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdir))) THEN
            got_var(idWdir)=.TRUE.
            S(ng)%Vid(idWdir)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdip))) THEN
            got_var(idWdip)=.TRUE.
            S(ng)%Vid(idWdip)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWpbt))) THEN
            got_var(idWpbt)=.TRUE.
            S(ng)%Vid(idWpbt)=var_id(i)
          END IF
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWorb).and.VarOut(idWorb,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWorb)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbrs).and.VarOut(idUbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.VarOut(idVbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.VarOut(idUbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.VarOut(idVbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.VarOut(idUbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.VarOut(idVbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVwc).and.VarOut(idUVwc,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUVwc)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.VarOut(idUbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.VarOut(idVbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.VarOut(idUbur,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.VarOut(idVbvr,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbvr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdir).and.VarOut(idWdir,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdir)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdip).and.VarOut(idWdip,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdip)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWpbt).and.VarOut(idWpbt,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWpbt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
      END IF QUERY
!
  10  FORMAT (/,' BBL_DEF_NF90 - unable to find variable: ',a,2x,       &
     &        ' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE bbl_def_nf90
!
!***********************************************************************
      SUBROUTINE bbl_wrt_nf90 (ng, model, tile,                         &
    &                          LBi, UBi, LBj, UBj,                      &
    &                          VarOut, S)
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
      integer :: gfactor, gtype, status
!
      real(dp) :: scale
!
      real(r8), allocatable :: wrk2d(:,:)
!
      character (len=256) :: ncname
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/BBL/bbl_output.F"//", bbl_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Bottom Boundary Layer model output variables into specified
!  output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 1737, MyFile)) RETURN
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
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (VarOut(idWorb,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWorb,             &
     &                       S(ng)%Vid(idWorb),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       FORCES(ng) % Uwave_rms)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWorb,             &
     &                       S(ng)%Vid(idWorb),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgWorb)
        END IF
        IF (FoundError(status, nf90_noerr, 1787, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWorb)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out current-induced, bottom U-stress at RHO-points.
!
      IF (VarOut(idUbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbrs,             &
     &                       S(ng)%Vid(idUbrs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % bustrc)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbrs,             &
     &                       S(ng)%Vid(idUbrs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgUbrs)
        END IF
        IF (FoundError(status, nf90_noerr, 1827, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbrs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out current-induced, bottom V-stress at RHO-points.
!
      IF (VarOut(idVbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbrs,             &
     &                       S(ng)%Vid(idVbrs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % bvstrc)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbrs,             &
     &                       S(ng)%Vid(idVbrs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgVbrs)
        END IF
        IF (FoundError(status, nf90_noerr, 1864, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbrs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom U-stress at RHO-points.
!
      IF (VarOut(idUbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbws,             &
     &                       S(ng)%Vid(idUbws),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % bustrw)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbws,             &
     &                       S(ng)%Vid(idUbws),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgUbws)
        END IF
        IF (FoundError(status, nf90_noerr, 1901, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom V-stress at RHO-points.
!
      IF (VarOut(idVbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbws,             &
     &                       S(ng)%Vid(idVbws),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % bvstrw)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbws,             &
     &                       S(ng)%Vid(idVbws),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgVbws)
        END IF
        IF (FoundError(status, nf90_noerr, 1938, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom U-stress at RHO-points.
!
      IF (VarOut(idUbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbcs,             &
     &                       S(ng)%Vid(idUbcs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % bustrcwmax)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbcs,             &
     &                       S(ng)%Vid(idUbcs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgUbcs)
        END IF
        IF (FoundError(status, nf90_noerr, 1975, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom V-stress at RHO-points.
!
      IF (VarOut(idVbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbcs,             &
     &                       S(ng)%Vid(idVbcs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % bvstrcwmax)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbcs,             &
     &                       S(ng)%Vid(idVbcs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgVbcs)
        END IF
        IF (FoundError(status, nf90_noerr, 2012, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wave and current bottom stress magnitude.
!
      IF (VarOut(idUVwc,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          IF (.not.allocated(wrk2d)) THEN
            allocate ( wrk2d(LBi:UBi, LBj:UBj) )
            wrk2d(LBi:UBi,LBj:UBj)=0.0_r8
          END IF
          wrk2d=SQRT(BBL(ng)%bustrcwmax*BBL(ng)%bustrcwmax+             &
     &               BBL(ng)%bvstrcwmax*BBL(ng)%bvstrcwmax+1.0E-10_r8)
!
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUVwc,             &
     &                       S(ng)%Vid(idUVwc),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       wrk2d)
          deallocate (wrk2d)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUVwc,             &
     &                       S(ng)%Vid(idUVwc),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgUVwc)
        END IF
        IF (FoundError(status, nf90_noerr, 2057, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVwc)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity at RHO-points.
!
      IF (VarOut(idUbot,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbot,             &
     &                       S(ng)%Vid(idUbot),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % Ubot)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbot,             &
     &                       S(ng)%Vid(idUbot),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgUbot)
        END IF
        IF (FoundError(status, nf90_noerr, 2094, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity at RHO-points
!
      IF (VarOut(idVbot,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbot,             &
     &                       S(ng)%Vid(idVbot),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % Vbot)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbot,             &
     &                       S(ng)%Vid(idVbot),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgVbot)
        END IF
        IF (FoundError(status, nf90_noerr, 2131, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-velocity above bed at RHO-points.
!
      IF (VarOut(idUbur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbur,             &
     &                       S(ng)%Vid(idUbur),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % Ur)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbur,             &
     &                       S(ng)%Vid(idUbur),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgUbur)
        END IF
        IF (FoundError(status, nf90_noerr, 2168, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbur)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-velocity above bed at RHO-points.
!
      IF (VarOut(idVbvr,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbvr,             &
     &                       S(ng)%Vid(idVbvr),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       BBL(ng) % Vr)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbvr,             &
     &                       S(ng)%Vid(idVbvr),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgVbvr)
        END IF
        IF (FoundError(status, nf90_noerr, 2204, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbvr)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced mean wave direction.
!
      IF (VarOut(idWdir,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdir,             &
     &                       S(ng)%Vid(idWdir),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       FORCES(ng) % Dwave)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdir,             &
     &                       S(ng)%Vid(idWdir),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgWdir)
        END IF
        IF (FoundError(status, nf90_noerr, 2441, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdir)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced peak wave direction.
!
      IF (VarOut(idWdip,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdip,             &
     &                       S(ng)%Vid(idWdip),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       FORCES(ng) % Dwavep)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdip,             &
     &                       S(ng)%Vid(idWdip),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgWdip)
        END IF
        IF (FoundError(status, nf90_noerr, 2480, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdip)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced bottom wave period.
!
      IF (VarOut(idWpbt,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWpbt,             &
     &                       S(ng)%Vid(idWpbt),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       FORCES(ng) % Pwave_bot)
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWpbt,             &
     &                       S(ng)%Vid(idWpbt),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgWpbt)
        END IF
        IF (FoundError(status, nf90_noerr, 2558, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWpbt)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
  10  FORMAT (/," BBL_WRT_NF90 - error while writing variable '",       &
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE bbl_wrt_nf90
!
      END MODULE bbl_output_mod
