/*
** kaarme.h - ROMS application options for external grid/boundary/river forcing
*/

/* Core dynamics */
#define SOLVE3D
#define SALINITY
#define SOLVE3D
#define AVERAGES

#define UV_ADV
#define UV_COR
#undef  UV_LDRAG
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define LIMIT_BSTRESS

#undef  TS_FIXED
#define TS_DIF2
#undef  TS_U3HADVECTION
#undef  TS_MPDATA
#define MIX_S_TS
#define DIAGNOSTICS
#undef  DJ_GRADPS
#undef  NONLIN_EOS
#undef  VAR_RHO_2D
#define DIAGNOSTICS_TS
#define DIAGNOSTICS_UV

#define MASKING
#undef  SALINITY

/* Disable nudging and sponge */
#undef  TCLM_NUDGING
#undef  M2CLM_NUDGING
#undef  M3CLM_NUDGING
#undef  ZCLM_NUDGING
#define SPONGE

/* Turbulence: use simpler BVF_MIXING initially for stability */
#undef  GLS_MIXING
#undef  KANTHA_CLAYSON
#undef  N2S2_HORAVG
#undef  RI_SPLINES
#undef  BVF_MIXING

/* Disable analytical grid; use external GRDNAME. Use analytic IC to avoid ININAME requirement */
#undef  ANA_GRID
#undef  ANA_INITIAL
#undef  ANA_SMFLUX
#undef  ANA_STFLUX
#undef  ANA_SSFLUX
#undef  ANA_BTFLUX
#undef  ANA_BSFLUX

/* Diagnostics can be enabled later if needed */
#undef  DIAGNOSTICS_TS
#undef  DIAGNOSTICS_UV
