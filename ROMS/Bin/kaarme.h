/*
** kaarme.h - ROMS application options for external grid/boundary/river forcing
*/

/* Core dynamics */
#define SOLVE3D
#define SALINITY
#define AVERAGES

#define UV_ADV
#define UV_COR
#define UV_LDRAG
#define UV_VIS2
#define MIX_S_UV

#define TS_DIF2
#define MIX_S_TS
#define DJ_GRADPS

#define MASKING

/* Turbulence: use GLS_MIXING (no analytical includes required) */
#define GLS_MIXING
#define KANTHA_CLAYSON
#define N2S2_HORAVG
#define RI_SPLINES

/* Disable analytical forcing/grid/IC; we use external files via roms.in */
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
