



















CBOP
C !ROUTINE: CPP_OPTIONS.h
C !INTERFACE:
C #include "CPP_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | main CPP options file for the model:
C | Control which optional features to compile in model/src code.
C *==================================================================*
CEOP

C CPP flags controlling particular source code features

C-- Forcing code options:

C o Shortwave heating as extra term in external_forcing.F
C Note: this should be a run-time option

C o Include/exclude Geothermal Heat Flux at the bottom of the ocean

C o Allow to account for heating due to friction (and momentum dissipation)

C o Allow mass source or sink of Fluid in the interior
C   (3-D generalisation of oceanic real-fresh water flux)

C o Include pressure loading code

C o Include/exclude balancing surface forcing fluxes code

C o Include/exclude balancing surface forcing relaxation code

C o Include/exclude checking for negative salinity

C-- Options to discard parts of the main code:

C o Exclude/allow external forcing-fields load
C   this allows to read & do simple linear time interpolation of oceanic
C   forcing fields, if no specific pkg (e.g., EXF) is used to compute them.

C o Include/exclude phi_hyd calculation code

C o Include/exclude sound speed calculation code
C o (Note that this is a diagnostic from Del Grasso algorithm, not derived from EOS)

C-- Vertical mixing code options:

C o Include/exclude call to S/R CONVECT

C o Include/exclude call to S/R CALC_DIFFUSIVITY

C o Allow full 3D specification of vertical diffusivity

C o Allow latitudinally varying BryanLewis79 vertical diffusivity

C o Exclude/allow partial-cell effect (physical or enhanced) in vertical mixing
C   this allows to account for partial-cell in vertical viscosity and diffusion,
C   either from grid-spacing reduction effect or as artificially enhanced mixing
C   near surface & bottom for too thin grid-cell

C o Exclude/allow to use isotropic 3-D Smagorinsky viscosity as diffusivity
C   for tracers (after scaling by constant Prandtl number)

C-- Time-stepping code options:

C o Include/exclude combined Surf.Pressure and Drag Implicit solver code

C o Include/exclude Implicit vertical advection code

C o Include/exclude AdamsBashforth-3rd-Order code

C o Include/exclude Quasi-Hydrostatic Stagger Time-step AdamsBashforth code

C-- Model formulation options:

C o Allow/exclude "Exact Convervation" of fluid in Free-Surface formulation
C   that ensures that d/dt(eta) is exactly equal to - Div.Transport

C o Allow the use of Non-Linear Free-Surface formulation
C   this implies that grid-cell thickness (hFactors) varies with time
C o Disable code for rStar coordinate and/or code for Sigma coordinate
c#define DISABLE_RSTAR_CODE
c#define DISABLE_SIGMA_CODE

C o Include/exclude nonHydrostatic code

C o Include/exclude GM-like eddy stress in momentum code

C-- Algorithm options:

C o Include/exclude code for Non Self-Adjoint (NSA) conjugate-gradient solver

C o Include/exclude code for single reduction Conjugate-Gradient solver

C o Choices for implicit solver routines solve_*diagonal.F
C   The following has low memory footprint, but not suitable for AD
C   The following one suitable for AD but does not vectorize

C-- Retired code options:

C o ALLOW isotropic scaling of harmonic and bi-harmonic terms when
C   using an locally isotropic spherical grid with (dlambda) x (dphi*cos(phi))
C *only for use on a lat-lon grid*
C   Setting this flag here affects both momentum and tracer equation unless
C   it is set/unset again in other header fields (e.g., GAD_OPTIONS.h).
C   The definition of the flag is commented to avoid interference with
C   such other header files.
C   The preferred method is specifying a value for viscAhGrid or viscA4Grid
C   in data which is then automatically scaled by the grid size;
C   the old method of specifying viscAh/viscA4 and this flag is provided
C   for completeness only (and for use with the adjoint).
c#define ISOTROPIC_COS_SCALING

C o This flag selects the form of COSINE(lat) scaling of bi-harmonic term.
C *only for use on a lat-lon grid*
C   Has no effect if ISOTROPIC_COS_SCALING is undefined.
C   Has no effect on vector invariant momentum equations.
C   Setting this flag here affects both momentum and tracer equation unless
C   it is set/unset again in other header fields (e.g., GAD_OPTIONS.h).
C   The definition of the flag is commented to avoid interference with
C   such other header files.
c#define COSINEMETH_III

C o Use "OLD" UV discretisation near boundaries (*not* recommended)
C   Note - only works with pkg/mom_fluxform and "no_slip_sides=.FALSE."
C          because the old code did not have no-slip BCs

C o Use LONG.bin, LATG.bin, etc., initialization for ini_curviliear_grid.F
C   Default is to use "new" grid files (OLD_GRID_IO undef) but OLD_GRID_IO
C   is still useful with, e.g., single-domain curvilinear configurations.

C o Use old EXTERNAL_FORCING_U,V,T,S subroutines (for backward compatibility)

C-- Other option files:

C o Execution environment support options
CBOP
C     !ROUTINE: CPP_EEOPTIONS.h
C     !INTERFACE:
C     include "CPP_EEOPTIONS.h"
C
C     !DESCRIPTION:
C     *==========================================================*
C     | CPP\_EEOPTIONS.h                                         |
C     *==========================================================*
C     | C preprocessor "execution environment" supporting        |
C     | flags. Use this file to set flags controlling the        |
C     | execution environment in which a model runs - as opposed |
C     | to the dynamical problem the model solves.               |
C     | Note: Many options are implemented with both compile time|
C     |       and run-time switches. This allows options to be   |
C     |       removed altogether, made optional at run-time or   |
C     |       to be permanently enabled. This convention helps   |
C     |       with the data-dependence analysis performed by the |
C     |       adjoint model compiler. This data dependency       |
C     |       analysis can be upset by runtime switches that it  |
C     |       is unable to recoginise as being fixed for the     |
C     |       duration of an integration.                        |
C     |       A reasonable way to use these flags is to          |
C     |       set all options as selectable at runtime but then  |
C     |       once an experimental configuration has been        |
C     |       identified, rebuild the code with the appropriate  |
C     |       options set at compile time.                       |
C     *==========================================================*
CEOP


C     In general the following convention applies:
C     ALLOW  - indicates an feature will be included but it may
C     CAN      have a run-time flag to allow it to be switched
C              on and off.
C              If ALLOW or CAN directives are "undef'd" this generally
C              means that the feature will not be available i.e. it
C              will not be included in the compiled code and so no
C              run-time option to use the feature will be available.
C
C     ALWAYS - indicates the choice will be fixed at compile time
C              so no run-time option will be present

C=== Macro related options ===
C--   Control storage of floating point operands
C     On many systems it improves performance only to use
C     8-byte precision for time stepped variables.
C     Constant in time terms ( geometric factors etc.. )
C     can use 4-byte precision, reducing memory utilisation and
C     boosting performance because of a smaller working set size.
C     However, on vector CRAY systems this degrades performance.
C     Enable to switch REAL4_IS_SLOW from genmake2 (with LET_RS_BE_REAL4):

C--   Control use of "double" precision constants.
C     Use D0 where it means REAL*8 but not where it means REAL*16

C=== IO related options ===
C--   Flag used to indicate whether Fortran formatted write
C     and read are threadsafe. On SGI the routines can be thread
C     safe, on Sun it is not possible - if you are unsure then
C     undef this option.

C--   Flag used to indicate whether Binary write to Local file (i.e.,
C     a different file for each tile) and read are thread-safe.

C--   Flag to turn off the writing of error message to ioUnit zero

C--   Alternative formulation of BYTESWAP, faster than
C     compiler flag -byteswapio on the Altix.

C--   Flag to turn on old default of opening scratch files with the
C     STATUS='SCRATCH' option. This method, while perfectly FORTRAN-standard,
C     caused filename conflicts on some multi-node/multi-processor platforms
C     in the past and has been replace by something (hopefully) more robust.

C--   Flag defined for eeboot_minimal.F, eeset_parms.F and open_copy_data_file.F
C     to write STDOUT, STDERR and scratch files from process 0 only.
C WARNING: to use only when absolutely confident that the setup is working
C     since any message (error/warning/print) from any proc <> 0 will be lost.

C=== MPI, EXCH and GLOBAL_SUM related options ===
C--   Flag turns off MPI_SEND ready_to_receive polling in the
C     gather_* subroutines to speed up integrations.

C--   Control MPI based parallel processing
CXXX We no longer select the use of MPI via this file (CPP_EEOPTIONS.h)
CXXX To use MPI, use an appropriate genmake2 options file or use
CXXX genmake2 -mpi .
CXXX #undef  ALLOW_USE_MPI

C--   Control use of communication that might overlap computation.
C     Under MPI selects/deselects "non-blocking" sends and receives.
C--   Control use of communication that is atomic to computation.
C     Under MPI selects/deselects "blocking" sends and receives.

C--   Control XY periodicity in processor to grid mappings
C     Note: Model code does not need to know whether a domain is
C           periodic because it has overlap regions for every box.
C           Model assume that these values have been
C           filled in some way.

C--   disconnect tiles (no exchange between tiles, just fill-in edges
C     assuming locally periodic subdomain)

C--   Always cumulate tile local-sum in the same order by applying MPI allreduce
C     to array of tiles ; can get slower with large number of tiles (big set-up)

C--   Alternative way of doing global sum without MPI allreduce call
C     but instead, explicit MPI send & recv calls. Expected to be slower.

C--   Alternative way of doing global sum on a single CPU
C     to eliminate tiling-dependent roundoff errors. Note: This is slow.

C=== Other options (to add/remove pieces of code) ===
C--   Flag to turn on checking for errors from all threads and procs
C     (calling S/R STOP_IF_ERROR) before stopping.

C--   Control use of communication with other component:
C     allow to import and export from/to Coupler interface.

C--   Activate some pieces of code for coupling to GEOS AGCM


CBOP
C     !ROUTINE: CPP_EEMACROS.h
C     !INTERFACE:
C     include "CPP_EEMACROS.h"
C     !DESCRIPTION:
C     *==========================================================*
C     | CPP_EEMACROS.h
C     *==========================================================*
C     | C preprocessor "execution environment" supporting
C     | macros. Use this file to define macros for  simplifying
C     | execution environment in which a model runs - as opposed
C     | to the dynamical problem the model solves.
C     *==========================================================*
CEOP


C     In general the following convention applies:
C     ALLOW  - indicates an feature will be included but it may
C     CAN      have a run-time flag to allow it to be switched
C              on and off.
C              If ALLOW or CAN directives are "undef'd" this generally
C              means that the feature will not be available i.e. it
C              will not be included in the compiled code and so no
C              run-time option to use the feature will be available.
C
C     ALWAYS - indicates the choice will be fixed at compile time
C              so no run-time option will be present

C     Flag used to indicate which flavour of multi-threading
C     compiler directives to use. Only set one of these.
C     USE_SOLARIS_THREADING  - Takes directives for SUN Workshop
C                              compiler.
C     USE_KAP_THREADING      - Takes directives for Kuck and
C                              Associates multi-threading compiler
C                              ( used on Digital platforms ).
C     USE_IRIX_THREADING     - Takes directives for SGI MIPS
C                              Pro Fortran compiler.
C     USE_EXEMPLAR_THREADING - Takes directives for HP SPP series
C                              compiler.
C     USE_C90_THREADING      - Takes directives for CRAY/SGI C90
C                              system F90 compiler.






C--   Define the mapping for the _BARRIER macro
C     On some systems low-level hardware support can be accessed through
C     compiler directives here.

C--   Define the mapping for the BEGIN_CRIT() and  END_CRIT() macros.
C     On some systems we simply execute this section only using the
C     master thread i.e. its not really a critical section. We can
C     do this because we do not use critical sections in any critical
C     sections of our code!

C--   Define the mapping for the BEGIN_MASTER_SECTION() and
C     END_MASTER_SECTION() macros. These are generally implemented by
C     simply choosing a particular thread to be "the master" and have
C     it alone execute the BEGIN_MASTER..., END_MASTER.. sections.

CcnhDebugStarts
C      Alternate form to the above macros that increments (decrements) a counter each
C      time a MASTER section is entered (exited). This counter can then be checked in barrier
C      to try and detect calls to BARRIER within single threaded sections.
C      Using these macros requires two changes to Makefile - these changes are written
C      below.
C      1 - add a filter to the CPP command to kill off commented _MASTER lines
C      2 - add a filter to the CPP output the converts the string N EWLINE to an actual newline.
C      The N EWLINE needs to be changes to have no space when this macro and Makefile changes
C      are used. Its in here with a space to stop it getting parsed by the CPP stage in these
C      comments.
C      #define IF ( a .EQ. 1 ) THEN  IF ( a .EQ. 1 ) THEN  N EWLINE      CALL BARRIER_MS(a)
C      #define ENDIF    CALL BARRIER_MU(a) N EWLINE        ENDIF
C      'CPP = cat $< | $(TOOLSDIR)/set64bitConst.sh |  grep -v '^[cC].*_MASTER' | cpp  -traditional -P'
C      .F.f:
C      $(CPP) $(DEFINES) $(INCLUDES) |  sed 's/N EWLINE/\n/' > $@
CcnhDebugEnds

C--   Control storage of floating point operands
C     On many systems it improves performance only to use
C     8-byte precision for time stepped variables.
C     Constant in time terms ( geometric factors etc.. )
C     can use 4-byte precision, reducing memory utilisation and
C     boosting performance because of a smaller working
C     set size. However, on vector CRAY systems this degrades
C     performance.
C- Note: global_sum/max macros were used to switch to  JAM routines (obsolete);
C  in addition, since only the R4 & R8 S/R are coded, GLOBAL RS & RL macros
C  enable to call the corresponding R4 or R8 S/R.



C- Note: a) exch macros were used to switch to  JAM routines (obsolete)
C        b) exch R4 & R8 macros are not practically used ; if needed,
C           will directly call the corrresponding S/R.

C--   Control use of JAM routines for Artic network (no longer supported)
C     These invoke optimized versions of "exchange" and "sum" that
C     utilize the programmable aspect of Artic cards.
CXXX No longer supported ; started to remove JAM routines.
CXXX #ifdef LETS_MAKE_JAM
CXXX #define CALL GLOBAL_SUM_R8 ( a, b) CALL GLOBAL_SUM_R8_JAM ( a, b)
CXXX #define CALL GLOBAL_SUM_R8 ( a, b ) CALL GLOBAL_SUM_R8_JAM ( a, b )
CXXX #define CALL EXCH_XY_RS ( a, b ) CALL EXCH_XY_R8_JAM ( a, b )
CXXX #define CALL EXCH_XY_RL ( a, b ) CALL EXCH_XY_R8_JAM ( a, b )
CXXX #define CALL EXCH_XYZ_RS ( a, b ) CALL EXCH_XYZ_R8_JAM ( a, b )
CXXX #define CALL EXCH_XYZ_RL ( a, b ) CALL EXCH_XYZ_R8_JAM ( a, b )
CXXX #endif

C--   Control use of "double" precision constants.
C     Use d0 where it means REAL*8 but not where it means REAL*16

C--   Substitue for 1.D variables
C     Sun compilers do not use 8-byte precision for literals
C     unless .Dnn is specified. CRAY vector machines use 16-byte
C     precision when they see .Dnn which runs very slowly!

C--   Set the format for writing processor IDs, e.g. in S/R eeset_parms
C     and S/R open_copy_data_file. The default of I9.9 should work for
C     a long time (until we will use 10e10 processors and more)



C o Include/exclude single header file containing multiple packages options
C   (AUTODIFF, COST, CTRL, ECCO, EXF ...) instead of the standard way where
C   each of the above pkg get its own options from its specific option file.
C   Although this method, inherited from ECCO setup, has been traditionally
C   used for all adjoint built, work is in progress to allow to use the
C   standard method also for adjoint built.
c#ifdef 
c# include "ECCO_CPPOPTIONS.h"
c#endif


C     Package-specific Options & Macros go here

C- use continuous power-law function for partition of energy between lateral
C  melting/freezing and thinning/thickening ; otherwise, use step function.

C- allow single grid-point debugging write to standard-output

C- only to check conservation
C  (change content of ICE_qleft,fresh,salFx-T files)

C CPP Macros go here


CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***

CBOP
C     !ROUTINE: THSICE_CALC_THICKN
C     !INTERFACE:
      SUBROUTINE THSICE_CALC_THICKN(
     I                  bi, bj,
     I                  iMin,iMax, jMin,jMax, dBugFlag,
     I                  iceMask, tFrz, tOce, v2oc,
     I                  snowP, prcAtm, sHeat, flxCnB,
     U                  icFrac, hIce, hSnow1, tSrf, qIc1, qIc2,
     U                  frwAtm, fzMlOc, flx2oc,
     O                  frw2oc, fsalt, frzSeaWat,
     I                  myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R  THSICE_CALC_THICKN
C     | o Calculate ice & snow thickness changes
C     *==========================================================*
C     \ev
C ADAPTED FROM:
C LANL CICE.v2.0.2
C-----------------------------------------------------------------------
C.. thermodynamics (vertical physics) based on M. Winton 3-layer model
C.. See Bitz, C. M. and W. H. Lipscomb, 1999:  An energy-conserving
C..       thermodynamic sea ice model for climate study.
C..       J. Geophys. Res., 104, 15669 - 15677.
C..     Winton, M., 1999:  "A reformulated three-layer sea ice model."
C..       Submitted to J. Atmos. Ocean. Technol.
C.. authors Elizabeth C. Hunke and William Lipscomb
C..         Fluid Dynamics Group, Los Alamos National Laboratory
C-----------------------------------------------------------------------
Cc****subroutine thermo_winton(n,fice,fsnow,dqice,dTsfc)
C.. Compute temperature change using Winton model with 2 ice layers, of
C.. which only the top layer has a variable heat capacity.

C---------------------------------
C  parameters that control the partitioning between lateral (ice area) and
C    vertical (ice thickness) ice volume changes.
C a) surface melting and bottom melting (runtime parameter: fracEnMelt):
C  frace is the fraction of available heat that is used for
C  lateral melting (and 1-frace reduces the thickness ) when
C o       hi < hThinIce        & frac > lowIcFrac2 : frace=1 (lateral melting only)
C o hThinIce < hi < hThickIce  & frac > lowIcFrac1 : frace=fracEnMelt
C o            hi > hThickIce or frac < lowIcFrac1 : frace=0 (thinning only)
C b) ocean freezing (and ice forming):
C - conductive heat flux (below sea-ice) always increases thickness.
C - under sea-ice, freezing potential (x iceFraction) is used to increase ice
C                  thickness or ice fraction (lateral growth), according to:
C o       hi < hThinIce       : use freezing potential to grow ice vertically;
C o hThinIce < hi < hThickIce : use partition factor fracEnFreez for lateral growth
c                               and (1-fracEnFreez) to increase thickness.
C o            hi > hThickIce : use all freezing potential to grow ice laterally
C                                (up to areaMax)
C - over open ocean, use freezing potential [x(1-iceFraction)] to grow ice laterally
C - lateral growth forms ice of the same or =hNewIceMax thickness, the less of the 2.
C - starts to form sea-ice over fraction iceMaskMin, as minimum ice-volume is reached
C---------------------------------
C     !USES:
      IMPLICIT NONE

C     == Global variables ===
C $Header: /u/gcmpack/MITgcm/verification/offline_exf_seaice/code/SIZE.h,v 1.4 2012/12/08 00:34:29 jmc Exp $
C $Name:  $

C
CBOP
C    !ROUTINE: SIZE.h
C    !INTERFACE:
C    include SIZE.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | SIZE.h Declare size of underlying computational grid.
C     *==========================================================*
C     | The design here support a three-dimensional model grid
C     | with indices I,J and K. The three-dimensional domain
C     | is comprised of nPx*nSx blocks of size sNx along one axis
C     | nPy*nSy blocks of size sNy along another axis and one
C     | block of size Nz along the final axis.
C     | Blocks have overlap regions of size OLx and OLy along the
C     | dimensions that are subdivided.
C     *==========================================================*
C     \ev
CEOP
C     Voodoo numbers controlling data layout.
C     sNx :: No. X points in sub-grid.
C     sNy :: No. Y points in sub-grid.
C     OLx :: Overlap extent in X.
C     OLy :: Overlat extent in Y.
C     nSx :: No. sub-grids in X.
C     nSy :: No. sub-grids in Y.
C     nPx :: No. of processes to use in X.
C     nPy :: No. of processes to use in Y.
C     Nx  :: No. points in X for the total domain.
C     Ny  :: No. points in Y for the total domain.
C     Nr  :: No. points in Z for full process domain.
      INTEGER sNx
      INTEGER sNy
      INTEGER OLx
      INTEGER OLy
      INTEGER nSx
      INTEGER nSy
      INTEGER nPx
      INTEGER nPy
      INTEGER Nx
      INTEGER Ny
      INTEGER Nr
      PARAMETER (
     &           sNx =  65,
     &           sNy =  65,
     &           OLx =   4,
     &           OLy =   4,
     &           nSx =   1,
     &           nSy =   1,
     &           nPx =   1,
     &           nPy =   1,
     &           Nx  = sNx*nSx*nPx,
     &           Ny  = sNy*nSy*nPy,
     &           Nr  =  1 )

C     MAX_OLX :: Set to the maximum overlap region size of any array
C     MAX_OLY    that will be exchanged. Controls the sizing of exch
C                routine buffers.
      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,
     &            MAX_OLY = OLy )

      INTEGER     nobcs
      PARAMETER ( nobcs = 4 )

CBOP
C     !ROUTINE: EEPARAMS.h
C     !INTERFACE:
C     include "EEPARAMS.h"
C
C     !DESCRIPTION:
C     *==========================================================*
C     | EEPARAMS.h                                               |
C     *==========================================================*
C     | Parameters for "execution environemnt". These are used   |
C     | by both the particular numerical model and the execution |
C     | environment support routines.                            |
C     *==========================================================*
CEOP

C     ========  EESIZE.h  ========================================

C     MAX_LEN_MBUF  :: Default message buffer max. size
C     MAX_LEN_FNAM  :: Default file name max. size
C     MAX_LEN_PREC  :: Default rec len for reading "parameter" files

      INTEGER MAX_LEN_MBUF
      PARAMETER ( MAX_LEN_MBUF = 512 )
      INTEGER MAX_LEN_FNAM
      PARAMETER ( MAX_LEN_FNAM = 512 )
      INTEGER MAX_LEN_PREC
      PARAMETER ( MAX_LEN_PREC = 200 )

C     MAX_NO_THREADS  :: Maximum number of threads allowed.
CC    MAX_NO_PROCS    :: Maximum number of processes allowed.
CC    MAX_NO_BARRIERS :: Maximum number of distinct thread "barriers"
      INTEGER MAX_NO_THREADS
      PARAMETER ( MAX_NO_THREADS =  4 )
c     INTEGER MAX_NO_PROCS
c     PARAMETER ( MAX_NO_PROCS   =  70000 )
c     INTEGER MAX_NO_BARRIERS
c     PARAMETER ( MAX_NO_BARRIERS = 1 )

C     Particularly weird and obscure voodoo numbers
C     lShare :: This wants to be the length in
C               [148]-byte words of the size of
C               the address "window" that is snooped
C               on an SMP bus. By separating elements in
C               the global sum buffer we can avoid generating
C               extraneous invalidate traffic between
C               processors. The length of this window is usually
C               a cache line i.e. small O(64 bytes).
C               The buffer arrays are usually short arrays
C               and are declared REAL ARRA(lShare[148],LBUFF).
C               Setting lShare[148] to 1 is like making these arrays
C               one dimensional.
      INTEGER cacheLineSize
      INTEGER lShare1
      INTEGER lShare4
      INTEGER lShare8
      PARAMETER ( cacheLineSize = 256 )
      PARAMETER ( lShare1 =  cacheLineSize )
      PARAMETER ( lShare4 =  cacheLineSize/4 )
      PARAMETER ( lShare8 =  cacheLineSize/8 )

CC    MAX_VGS  :: Maximum buffer size for Global Vector Sum
c     INTEGER MAX_VGS
c     PARAMETER ( MAX_VGS = 8192 )

C     ========  EESIZE.h  ========================================

C     Symbolic values
C     precXXXX :: precision used for I/O
      INTEGER precFloat32
      PARAMETER ( precFloat32 = 32 )
      INTEGER precFloat64
      PARAMETER ( precFloat64 = 64 )

C     Real-type constant for some frequently used simple number (0,1,2,1/2):
      Real*8     zeroRS, oneRS, twoRS, halfRS
      PARAMETER ( zeroRS = 0.0D0 , oneRS  = 1.0D0 )
      PARAMETER ( twoRS  = 2.0D0 , halfRS = 0.5D0 )
      Real*8     zeroRL, oneRL, twoRL, halfRL
      PARAMETER ( zeroRL = 0.0D0 , oneRL  = 1.0D0 )
      PARAMETER ( twoRL  = 2.0D0 , halfRL = 0.5D0 )

C     UNSET_xxx :: Used to indicate variables that have not been given a value
      Real*8  UNSET_FLOAT8
      PARAMETER ( UNSET_FLOAT8 = 1.234567D5 )
      Real*4  UNSET_FLOAT4
      PARAMETER ( UNSET_FLOAT4 = 1.234567E5 )
      Real*8     UNSET_RL
      PARAMETER ( UNSET_RL     = 1.234567D5 )
      Real*8     UNSET_RS
      PARAMETER ( UNSET_RS     = 1.234567D5 )
      INTEGER UNSET_I
      PARAMETER ( UNSET_I      = 123456789  )

C     debLevX  :: used to decide when to print debug messages
      INTEGER debLevZero
      INTEGER debLevA, debLevB,  debLevC, debLevD, debLevE
      PARAMETER ( debLevZero=0 )
      PARAMETER ( debLevA=1 )
      PARAMETER ( debLevB=2 )
      PARAMETER ( debLevC=3 )
      PARAMETER ( debLevD=4 )
      PARAMETER ( debLevE=5 )

C     SQUEEZE_RIGHT      :: Flag indicating right blank space removal
C                           from text field.
C     SQUEEZE_LEFT       :: Flag indicating left blank space removal
C                           from text field.
C     SQUEEZE_BOTH       :: Flag indicating left and right blank
C                           space removal from text field.
C     PRINT_MAP_XY       :: Flag indicating to plot map as XY slices
C     PRINT_MAP_XZ       :: Flag indicating to plot map as XZ slices
C     PRINT_MAP_YZ       :: Flag indicating to plot map as YZ slices
C     commentCharacter   :: Variable used in column 1 of parameter
C                           files to indicate comments.
C     INDEX_I            :: Variable used to select an index label
C     INDEX_J               for formatted input parameters.
C     INDEX_K
C     INDEX_NONE
      CHARACTER*(*) SQUEEZE_RIGHT
      PARAMETER ( SQUEEZE_RIGHT = 'R' )
      CHARACTER*(*) SQUEEZE_LEFT
      PARAMETER ( SQUEEZE_LEFT = 'L' )
      CHARACTER*(*) SQUEEZE_BOTH
      PARAMETER ( SQUEEZE_BOTH = 'B' )
      CHARACTER*(*) PRINT_MAP_XY
      PARAMETER ( PRINT_MAP_XY = 'XY' )
      CHARACTER*(*) PRINT_MAP_XZ
      PARAMETER ( PRINT_MAP_XZ = 'XZ' )
      CHARACTER*(*) PRINT_MAP_YZ
      PARAMETER ( PRINT_MAP_YZ = 'YZ' )
      CHARACTER*(*) commentCharacter
      PARAMETER ( commentCharacter = '#' )
      INTEGER INDEX_I
      INTEGER INDEX_J
      INTEGER INDEX_K
      INTEGER INDEX_NONE
      PARAMETER ( INDEX_I    = 1,
     &            INDEX_J    = 2,
     &            INDEX_K    = 3,
     &            INDEX_NONE = 4 )

C     EXCH_IGNORE_CORNERS :: Flag to select ignoring or
C     EXCH_UPDATE_CORNERS    updating of corners during an edge exchange.
      INTEGER EXCH_IGNORE_CORNERS
      INTEGER EXCH_UPDATE_CORNERS
      PARAMETER ( EXCH_IGNORE_CORNERS = 0,
     &            EXCH_UPDATE_CORNERS = 1 )

C     FORWARD_SIMULATION
C     REVERSE_SIMULATION
C     TANGENT_SIMULATION
      INTEGER FORWARD_SIMULATION
      INTEGER REVERSE_SIMULATION
      INTEGER TANGENT_SIMULATION
      PARAMETER ( FORWARD_SIMULATION = 0,
     &            REVERSE_SIMULATION = 1,
     &            TANGENT_SIMULATION = 2 )

C--   COMMON /EEPARAMS_L/ Execution environment public logical variables.
C     eeBootError    :: Flags indicating error during multi-processing
C     eeEndError     :: initialisation and termination.
C     fatalError     :: Flag used to indicate that the model is ended with an error
C     debugMode      :: controls printing of debug msg (sequence of S/R calls).
C     useSingleCpuIO :: When useSingleCpuIO is set, MDS_WRITE_FIELD outputs from
C                       master MPI process only. -- NOTE: read from main parameter
C                       file "data" and not set until call to INI_PARMS.
C     useSingleCpuInput :: When useSingleCpuInput is set, EXF_INTERP_READ
C                       reads forcing files from master MPI process only.
C                       -- NOTE: read from main parameter file "data"
C                          and defaults to useSingleCpuInput = useSingleCpuIO
C     printMapIncludesZeros  :: Flag that controls whether character constant
C                               map code ignores exact zero values.
C     useCubedSphereExchange :: use Cubed-Sphere topology domain.
C     useCoupler     :: use Coupler for a multi-components set-up.
C     useNEST_PARENT :: use Parent Nesting interface (pkg/nest_parent)
C     useNEST_CHILD  :: use Child  Nesting interface (pkg/nest_child)
C     useNest2W_parent :: use Parent 2-W Nesting interface (pkg/nest2w_parent)
C     useNest2W_child  :: use Child  2-W Nesting interface (pkg/nest2w_child)
C     useOASIS       :: use OASIS-coupler for a multi-components set-up.
      COMMON /EEPARAMS_L/
c    &  eeBootError, fatalError, eeEndError,
     &  eeBootError, eeEndError, fatalError, debugMode,
     &  useSingleCpuIO, useSingleCpuInput, printMapIncludesZeros,
     &  useCubedSphereExchange, useCoupler,
     &  useNEST_PARENT, useNEST_CHILD,
     &  useNest2W_parent, useNest2W_child, useOASIS,
     &  useSETRLSTK, useSIGREG
      LOGICAL eeBootError
      LOGICAL eeEndError
      LOGICAL fatalError
      LOGICAL debugMode
      LOGICAL useSingleCpuIO
      LOGICAL useSingleCpuInput
      LOGICAL printMapIncludesZeros
      LOGICAL useCubedSphereExchange
      LOGICAL useCoupler
      LOGICAL useNEST_PARENT
      LOGICAL useNEST_CHILD
      LOGICAL useNest2W_parent
      LOGICAL useNest2W_child
      LOGICAL useOASIS
      LOGICAL useSETRLSTK
      LOGICAL useSIGREG

C--   COMMON /EPARAMS_I/ Execution environment public integer variables.
C     errorMessageUnit    :: Fortran IO unit for error messages
C     standardMessageUnit :: Fortran IO unit for informational messages
C     maxLengthPrt1D :: maximum length for printing (to Std-Msg-Unit) 1-D array
C     scrUnit1      :: Scratch file 1 unit number
C     scrUnit2      :: Scratch file 2 unit number
C     eeDataUnit    :: Unit # for reading "execution environment" parameter file
C     modelDataUnit :: Unit number for reading "model" parameter file.
C     numberOfProcs :: Number of processes computing in parallel
C     pidIO         :: Id of process to use for I/O.
C     myBxLo, myBxHi :: Extents of domain in blocks in X and Y
C     myByLo, myByHi :: that each threads is responsble for.
C     myProcId      :: My own "process" id.
C     myPx          :: My X coord on the proc. grid.
C     myPy          :: My Y coord on the proc. grid.
C     myXGlobalLo   :: My bottom-left (south-west) x-index global domain.
C                      The x-coordinate of this point in for example m or
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C     myYGlobalLo   :: My bottom-left (south-west) y-index in global domain.
C                      The y-coordinate of this point in for example m or
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C     nThreads      :: No. of threads
C     nTx, nTy      :: No. of threads in X and in Y
C                      This assumes a simple cartesian gridding of the threads
C                      which is not required elsewhere but that makes it easier
C     ioErrorCount  :: IO Error Counter. Set to zero initially and increased
C                      by one every time an IO error occurs.
      COMMON /EEPARAMS_I/
     &  errorMessageUnit, standardMessageUnit, maxLengthPrt1D,
     &  scrUnit1, scrUnit2, eeDataUnit, modelDataUnit,
     &  numberOfProcs, pidIO, myProcId,
     &  myPx, myPy, myXGlobalLo, myYGlobalLo, nThreads,
     &  myBxLo, myBxHi, myByLo, myByHi,
     &  nTx, nTy, ioErrorCount
      INTEGER errorMessageUnit
      INTEGER standardMessageUnit
      INTEGER maxLengthPrt1D
      INTEGER scrUnit1
      INTEGER scrUnit2
      INTEGER eeDataUnit
      INTEGER modelDataUnit
      INTEGER ioErrorCount(MAX_NO_THREADS)
      INTEGER myBxLo(MAX_NO_THREADS)
      INTEGER myBxHi(MAX_NO_THREADS)
      INTEGER myByLo(MAX_NO_THREADS)
      INTEGER myByHi(MAX_NO_THREADS)
      INTEGER myProcId
      INTEGER myPx
      INTEGER myPy
      INTEGER myXGlobalLo
      INTEGER myYGlobalLo
      INTEGER nThreads
      INTEGER nTx
      INTEGER nTy
      INTEGER numberOfProcs
      INTEGER pidIO

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
C     *==========================================================*
C     | THSICE_SIZE.h Declare size of arrays for Therm_SeaIce pkg
C     *==========================================================*

C.. number layers of ice
C     nlyr   ::   maximum number of ice layers
      INTEGER nlyr
      PARAMETER (nlyr = 2)

C--   Energy distribution (lateral / thickening-thinning) using a power law:
C     power-law exponent is set to: 1+2^powerLawExp2
      INTEGER     powerLawExp2
      PARAMETER ( powerLawExp2 = 2 )

C--   identifiers for advected properties
      INTEGER GAD_SI_FRAC, GAD_SI_HSNOW
      INTEGER GAD_SI_HICE, GAD_SI_QICE1, GAD_SI_QICE2
      PARAMETER ( GAD_SI_FRAC  = -5,
     &            GAD_SI_HSNOW = -6,
     &            GAD_SI_HICE  = -7,
     &            GAD_SI_QICE1 = -8,
     &            GAD_SI_QICE2 = -9 )



CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
C     *==========================================================*
C     | THSICE_PARAMS.h
C     | o Header file for Therm_SeaIce package parameters:
C     |   - basic parameter ( I/O frequency, etc ...)
C     |   - physical constants (used in therm_SeaIce pkg)
C     *==========================================================*

C----------------------------------------------------------------------------
C.. Common blocks for almost everything that the sea ice model passes around.
C----------------------------------------------------------------------------

C--   COMMON / THSICE_PHYSPAR_R / physical (real) parameter
C.. densities
C     rhos      ::   density of snow [kg/m^3]
C     rhoi      ::   density of ice [kg/m^3]
C     rhosw     ::   density of seawater [kg/m^3]
C     rhofw     ::   density of fresh water [kg/m^3]
C     floodFac  ::   flooding factor = (rhosw-rhoi)/rhos [dimensionless]
C.. specific heats
C     cpIce     ::   specific heat of fresh ice [J/kg/K]
C     cpWater   ::   specific heat of water [J/kg/K]
C .. thermal conductivity. QQ check units
C     kIce      ::   thermal conductivity of pure ice [W/m/K]
C     kSnow     ::   thermal conductivity of snow [W/m/K]
C .. heat transfer coefficient
C     bMeltCoef ::   base-melting heat transfer coefficient
C                    (between ice & water) [no unit]
C .. latent heat
C     Lfresh    ::   latent heat of melting of pure ice [J/kg]
C .. Enthalpy
C     qsnow     ::   snow enthalpy [J/kg]
C .. Albedo
C     albColdSnow :: albedo of cold (=dry) new snow (Tsfc < tempSnowAlb)
C     albWarmSnow :: albedo of warm (=wet) new snow (Tsfc = 0)
C     tempSnowAlb :: temperature transition from ColdSnow to WarmSnow Alb. [oC]
C     albOldSnow  :: albedo of old snow (snowAge > 35.d)
C     albIceMax   :: max albedo of bare ice (thick ice)
C     albIceMin   :: minimum ice albedo (very thin ice)
C     hAlbIce     :: ice thickness for albedo transition: thin/thick ice albedo
C     hAlbSnow    :: snow thickness for albedo transition: snow/ice albedo
C     hNewSnowAge :: new snow thickness that refresh the snow-age (by 1/e)
C     snowAgTime  :: snow aging time scale [s]
C .. Solar parameters
C     i0swFrac  ::   fraction of penetrating solar rad
C     ksolar    ::   bulk solar abs coeff of sea ice [m^-1]
C     dhSnowLin ::   half slope of linear distribution of snow thickness within
C                    the grid-cell (from hSnow-dhSnow to hSnow+dhSnow, if full
C                    ice & snow cover) [m] ; (only used for SW radiation).
C .. Salinity
C     saltIce   ::   salinity of ice [g/kg]
C     S_winton  ::   Winton salinity of ice [g/kg]
C .. melting
C     Tf0kel    ::   Freezing temp of fresh water in Kelvin = 273.15
C     mu_Tf     ::   linear dependence of melting temperature on salinity [oC/(g/kg)]
C                     Tf(sea-water) = -mu_Tf * S
C     Tmlt1     ::   Winton melting temperature: Tmlt1 = -mu_Tf * S_winton
C     Terrmax   ::   Temperature convergence criteria [oC]
C .. Min/Max
C     hIceMin   ::   Minimum ice  thickness [m]
C     hiMax     ::   Maximum ice  thickness [m]
C     hsMax     ::   Maximum snow thickness [m]
C .. for fractional ice
C     iceMaskMax  :: maximum Ice fraction (=1 for no fractional ice)
C     iceMaskMin  :: mimimum Ice fraction (=1 for no fractional ice)
C     fracEnFreez :: fraction of energy going to lateral freezing (vs height increase)
C     fracEnMelt  :: fraction of energy going to lateral melting (vs height decrease)
C                                         (=0 for no fract. ice)
C     hThinIce    :: ice height above which fracEnMelt/Freez are applied [m]
C                                         (=hIceMin for no fractional ice)
C     hThickIce   :: ice height below which fracEnMelt/Freez are applied [m]
C                                         (=large for no fractional ice)
C     hNewIceMax  :: new ice maximum thickness [m]
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      COMMON / THSICE_PHYSPAR_R /
     &  rhos, rhoi, rhosw, rhofw, floodFac,
     &  cpIce, cpWater,
     &  kIce, kSnow,
     &  bMeltCoef, Lfresh, qsnow,
     &  albColdSnow, albWarmSnow, tempSnowAlb,
     &  albOldSnow, hNewSnowAge, snowAgTime,
     &  albIceMax, albIceMin, hAlbIce, hAlbSnow,
     &  i0swFrac, ksolar, dhSnowLin,
     &  saltIce, S_winton, mu_Tf,
     &  Tf0kel, Tmlt1, Terrmax,
     &  hIceMin, hiMax, hsMax,
     &  iceMaskMax, iceMaskMin,
     &  fracEnMelt, fracEnFreez,
     &  hThinIce, hThickIce, hNewIceMax

      Real*8  rhos
      Real*8  rhoi
      Real*8  rhosw
      Real*8  rhofw
      Real*8  floodFac
      Real*8  cpIce
      Real*8  cpWater
      Real*8  kIce
      Real*8  kSnow
      Real*8  bMeltCoef
      Real*8  Lfresh
      Real*8  qsnow
      Real*8  albColdSnow
      Real*8  albWarmSnow
      Real*8  tempSnowAlb
      Real*8  albOldSnow
      Real*8  hNewSnowAge
      Real*8  snowAgTime
      Real*8  albIceMax
      Real*8  albIceMin
      Real*8  hAlbIce
      Real*8  hAlbSnow
      Real*8  i0swFrac
      Real*8  ksolar
      Real*8  dhSnowLin
      Real*8  saltIce
      Real*8  S_winton
      Real*8  mu_Tf
      Real*8  Tf0kel
      Real*8  Tmlt1
      Real*8  Terrmax
      Real*8  hIceMin
      Real*8  hiMax
      Real*8  hsMax
      Real*8 iceMaskMax
      Real*8 iceMaskMin
      Real*8 fracEnMelt
      Real*8 fracEnFreez
      Real*8 hThinIce
      Real*8 hThickIce
      Real*8 hNewIceMax

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   COMMON / THSICE_PAR_L / ice model (logical) parameters
C     stepFwd_oceMxL        :: step forward mixed-layer T & S (slab-ocean)
C     thSIce_skipThermo     :: by-pass seaice thermodynamics
C     thSIce_calc_albNIR    :: calculate Near Infra-Red Albedo
C     thSIce_tave_mdsio     :: write TimeAverage output using MDSIO
C     thSIce_snapshot_mdsio :: write snap-shot output   using MDSIO
C     thSIce_mon_stdio      :: write monitor to std-outp
C     thSIce_tave_mnc       :: write TimeAverage output using MNC
C     thSIce_snapshot_mnc   :: write snap-shot output   using MNC
C     thSIce_mon_mnc        :: write monitor to netcdf file
C     thSIce_pickup_read_mnc    :: pickup read w/ MNC
C     thSIce_pickup_write_mnc   :: pickup write w/ MNC
C     thSIce_pickup_write_mdsio :: pickup write w/ MDSIO
      COMMON / THSICE_PAR_L /
     &     stepFwd_oceMxL, thSIce_skipThermo,
     &     thSIce_calc_albNIR,
     &     thSIce_tave_mdsio, thSIce_snapshot_mdsio, thSIce_mon_stdio,
     &     thSIce_tave_mnc,   thSIce_snapshot_mnc,   thSIce_mon_mnc,
     &     thSIce_pickup_read_mnc,
     &     thSIce_pickup_write_mdsio,
     &     thSIce_pickup_write_mnc

      LOGICAL stepFwd_oceMxL
      LOGICAL thSIce_skipThermo
      LOGICAL thSIce_calc_albNIR
      LOGICAL thSIce_tave_mdsio, thSIce_snapshot_mdsio, thSIce_mon_stdio
      LOGICAL thSIce_tave_mnc,   thSIce_snapshot_mnc,   thSIce_mon_mnc
      LOGICAL thSIce_pickup_read_mnc
      LOGICAL thSIce_pickup_write_mdsio
      LOGICAL thSIce_pickup_write_mnc

C--   COMMON / THSICE_PAR_I / ice model (integer) parameters
C     startIceModel   :: =1 : start ice model at nIter0 ; =0 : use pickup files
C                     :: -1 : start from a small pickup (without Mix.Layer)
C     nitMaxTsf       :: maximum Nb of iter to find Surface Temp (Trsf)
C     thSIceAdvScheme :: thSIce Advection scheme selector
C     thSIceBalanceAtmFW :: select balancing Fresh-Water flux from Atm+Land
C                        :: =0 : none ; =1 : uniform ; =2 : scaled by Precip
      COMMON / THSICE_PAR_I /
     &  startIceModel, nitMaxTsf, thSIceAdvScheme, thSIceBalanceAtmFW

      INTEGER startIceModel
      INTEGER nitMaxTsf
      INTEGER thSIceAdvScheme
      INTEGER thSIceBalanceAtmFW

C--   COMMON / THSICE_PAR_R / ice model (real) parameters
C     thSIce_deltaT   :: ice model time-step, seaice thicken/extend [s]
C     thSIce_dtTemp   :: ice model time-step, solve4temp [s]
C     ocean_deltaT    :: ocean mixed-layer time-step [s]
C     tauRelax_MxL    :: Relaxation time scale for MixLayer T [s]
C     tauRelax_MxL_salt :: Relaxation time scale for MixLayer S [s]
C     hMxL_default    :: default value for ocean MixLayer thickness [m]
C     sMxL_default    :: default value for salinity in MixLayer [g/kg]
C     vMxL_default    :: default value for ocean current velocity in MxL [m/s]
C     thSIce_diffK    :: thickness (horizontal) diffusivity [m^2/s]
C     stressReduction :: reduction factor for wind-stress under sea-ice [0-1]
C     thSIce_taveFreq :: Frequency^-1 for time-Aver. output [s]
C     thSIce_diagFreq :: Frequency^-1 for diagnostic output [s]
C     thSIce_monFreq  :: Frequency^-1 for monitor    output [s]
      COMMON / THSICE_PAR_R /
     &  thSIce_deltaT,  thSIce_dtTemp, ocean_deltaT,
     &  tauRelax_MxL, tauRelax_MxL_salt,
     &  hMxL_default,  sMxL_default, vMxL_default,
     &  thSIce_diffK,  stressReduction,
     &  thSIce_taveFreq, thSIce_diagFreq, thSIce_monFreq

      Real*8  thSIce_deltaT, thSIce_dtTemp, ocean_deltaT
      Real*8  tauRelax_MxL, tauRelax_MxL_salt
      Real*8  hMxL_default, sMxL_default, vMxL_default
      Real*8  thSIce_diffK,  stressReduction
      Real*8  thSIce_taveFreq, thSIce_diagFreq, thSIce_monFreq

C--   COMMON / THSICE_PAR_C / ice model (character) parameters
C     thSIceFract_InitFile :: File name for initial ice fraction
C     thSIceThick_InitFile :: File name for initial ice thickness
C     thSIceSnowH_InitFile :: File name for initial snow thickness
C     thSIceSnowA_InitFile :: File name for initial snow Age
C     thSIceEnthp_InitFile :: File name for initial ice enthalpy
C     thSIceTsurf_InitFile :: File name for initial surf. temp
      COMMON / THSICE_PAR_C /
     &  thSIceFract_InitFile,
     &  thSIceThick_InitFile,
     &  thSIceSnowH_InitFile,
     &  thSIceSnowA_InitFile,
     &  thSIceEnthp_InitFile,
     &  thSIceTsurf_InitFile
      CHARACTER*(MAX_LEN_FNAM) thSIceFract_InitFile
      CHARACTER*(MAX_LEN_FNAM) thSIceThick_InitFile
      CHARACTER*(MAX_LEN_FNAM) thSIceSnowH_InitFile
      CHARACTER*(MAX_LEN_FNAM) thSIceSnowA_InitFile
      CHARACTER*(MAX_LEN_FNAM) thSIceEnthp_InitFile
      CHARACTER*(MAX_LEN_FNAM) thSIceTsurf_InitFile

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine Arguments ==
C     bi,bj       :: tile indices
C     iMin,iMax   :: computation domain: 1rst index range
C     jMin,jMax   :: computation domain: 2nd  index range
C     dBugFlag    :: allow to print debugging stuff (e.g. on 1 grid point).
C---  Input:
C         iceMask :: sea-ice fractional mask [0-1]
C  tFrz           :: sea-water freezing temperature [oC] (function of S)
C  tOce           :: surface level oceanic temperature [oC]
C  v2oc           :: square of ocean surface-level velocity [m2/s2]
C  snowP          :: snow precipitation                [kg/m2/s]
C  prcAtm         :: total precip from the atmosphere [kg/m2/s]
C  sHeat          :: surf heating flux left to melt snow or ice (= Atmos-conduction)
C  flxCnB         :: heat flux conducted through the ice to bottom surface
C---  Modified (input&output):
C  icFrac         :: fraction of grid area covered in ice
C  hIce           :: ice height [m]
C  hSnow1          :: snow height [m]
C  tSrf           :: surface (ice or snow) temperature
C  qIc1   (qicen) :: ice enthalpy (J/kg), 1rst level
C  qIc2   (qicen) :: ice enthalpy (J/kg), 2nd level
C  frwAtm (evpAtm):: evaporation to the atmosphere [kg/m2/s] (>0 if evaporate)
C  fzMlOc         :: ocean mixed-layer freezing/melting potential [W/m2]
C  flx2oc         :: net heat flux to ocean    [W/m2]          (> 0 downward)
C---  Output
C  frw2oc         :: Total fresh water flux to ocean [kg/m2/s] (> 0 downward)
C  fsalt          :: salt flux to ocean        [g/m2/s]        (> 0 downward)
C  frzSeaWat      :: seawater freezing rate (expressed as mass flux) [kg/m^2/s]
C---  Input:
C     myTime      :: current Time of simulation [s]
C     myIter      :: current Iteration number in simulation
C     myThid      :: my Thread Id number
      INTEGER bi,bj
      INTEGER iMin, iMax
      INTEGER jMin, jMax
      LOGICAL dBugFlag
      Real*8 iceMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 tFrz   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 tOce   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 v2oc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 snowP  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 prcAtm (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 sHeat  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 flxCnB (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 icFrac (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 hIce   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 hSnow1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 tSrf   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 qIc1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 qIc2   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 frwAtm (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 fzMlOc (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 flx2oc (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 frw2oc (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 fsalt  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 frzSeaWat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  myTime
      INTEGER myIter
      INTEGER myThid
CEOP


C     !LOCAL VARIABLES:
C---  local copy of input/output argument list variables (see description above)
      Real*8 qicen(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nlyr)

C     == Local Variables ==
C     i,j,k      :: loop indices
C     rec_nlyr   :: reciprocal of number of ice layers (real value)
C     evapLoc    :: evaporation over snow/ice [kg/m2/s] (>0 if evaporate)
C     Fbot       :: oceanic heat flux used to melt/form ice [W/m2]
C     etop       :: energy for top melting    (J m-2)
C     ebot       :: energy for bottom melting (J m-2)
C     etope      :: energy (from top)    for lateral melting (J m-2)
C     ebote      :: energy (from bottom) for lateral melting (J m-2)
C     extend     :: total energy for lateral melting (J m-2)
C     hnew(nlyr) :: new ice layer thickness (m)
C     hlyr       :: individual ice layer thickness (m)
C     dhi        :: change in ice thickness
C     dhs        :: change in snow thickness
C     rq         :: rho * q for a layer
C     rqh        :: rho * q * h for a layer
C     qbot       :: enthalpy for new ice at bottom surf (J/kg)
C     dt         :: timestep
C     esurp      :: surplus energy from melting (J m-2)
C     mwater0    :: fresh water mass gained/lost (kg/m^2)
C     msalt0     :: salt gained/lost  (kg/m^2)
C     freshe     :: fresh water gain from extension melting
C     salte      :: salt gained from extension melting
C     lowIcFrac1 :: ice-fraction lower limit above which partial (lowIcFrac1)
C     lowIcFrac2 :: or full (lowIcFrac2) lateral melting is allowed.
C     from THSICE_RESHAPE_LAYERS
C     f1         :: Fraction of upper layer ice in new layer
C     qh1, qh2   :: qice*h for layers 1 and 2
C     qhtot      :: qh1 + qh2
C     q2tmp      :: Temporary value of qice for layer 2
      INTEGER  i,j,k
      Real*8 rec_nlyr
      Real*8 evapLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 Fbot   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 etop   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 ebot   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 etope  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 ebote  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 esurp  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 extend
      Real*8 hnew   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nlyr)
      Real*8 hlyr
      Real*8 dhi
      Real*8 dhs
      Real*8 rq
      Real*8 rqh
      Real*8 qbot
      Real*8 dt
      Real*8 mwater0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 msalt0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 freshe
      Real*8 salte
      Real*8 lowIcFrac1, lowIcFrac2
      Real*8  f1
      Real*8  qh1, qh2
      Real*8  qhtot
      Real*8  q2tmp

      Real*8  ustar, cpchr
      Real*8  chi
      Real*8  frace, rs, hq
      INTEGER powerLaw
      Real*8 rec_pLaw
      Real*8 c1Mlt, c2Mlt, aMlt, hMlt
      Real*8 c1Frz, c2Frz, aFrz, hFrz
      Real*8 enFrcMlt(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 xxMlt, tmpMlt
      Real*8 enFrcFrz(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 xxFrz, tmpFrz

C-    define grid-point location where to print debugging values
CBOP
C !ROUTINE: THSICE_DEBUG.h

C !INTERFACE:
C #include "THSICE_DEBUG.h"
C       LOGICAL dBug
C       dBug(i,j,bi,bj)

C !DESCRIPTION:
C Function used for debugging: define a single grid-point location
C  where values of various variables are printed (to standard-output file)
C
CEOP

      LOGICAL dBug
c     dBug(i,j,bi,bj) = .FALSE.
      dBug(i,j,bi,bj) = dBugFlag .AND.
     &         ( i.EQ.15 .AND. j.EQ.11 .AND. bi.EQ.3 .AND. bj.EQ.1 )

 1020 FORMAT(A,1P4E11.3)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


      rec_nlyr = nlyr
      rec_nlyr = 1.D0 / rec_nlyr
      dt  = thSIce_deltaT

C     for now, use hard coded threshold (iceMaskMin +1.% and +10.%)
      lowIcFrac1 = iceMaskMin*1.01D0
      lowIcFrac2 = iceMaskMin*1.10D0
      IF ( powerLawExp2 .GE. 1 ) THEN
        powerLaw = 1 + 2**powerLawExp2
        rec_pLaw = powerLaw
        rec_pLaw = 1.D0 / rec_pLaw
C-    Coef for melting:
C     lateral-melting energy fraction = fracEnMelt - [ aMlt*(hi-hMlt) ]^powerLaw
        c1Mlt = fracEnMelt**rec_pLaw
        c2Mlt = (1.D0 - fracEnMelt)**rec_pLaw
        aMlt = (c1Mlt+c2Mlt)/(hThickIce-hThinIce)
        hMlt = hThinIce+c2Mlt/aMlt
C-    Coef for freezing:
C     thickening energy fraction     = fracEnFreez - [ aFrz*(hi-hFrz) ]^powerLaw
        c1Frz = fracEnFreez**rec_pLaw
        c2Frz = (1.D0 -fracEnFreez)**rec_pLaw
        aFrz = (c1Frz+c2Frz)/(hThickIce-hThinIce)
        hFrz = hThinIce+c2Frz/aFrz
      ELSE
C-    Linear relation
        powerLaw = 1
        aMlt = -1.D0 /(hThickIce-hThinIce)
        hMlt = hThickIce
        aFrz = -1.D0 /(hThickIce-hThinIce)
        hFrz = hThickIce
      ENDIF

C     initialise local arrays
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        evapLoc(i,j) = 0.D0
        Fbot   (i,j) = 0.D0
        etop   (i,j) = 0.D0
        ebot   (i,j) = 0.D0
        etope  (i,j) = 0.D0
        ebote  (i,j) = 0.D0
        esurp  (i,j) = 0.D0
        mwater0(i,j) = 0.D0
        msalt0 (i,j) = 0.D0
        enFrcMlt(i,j)= 0.D0
        enFrcFrz(i,j)= 0.D0
       ENDDO
      ENDDO
      DO k = 1,nlyr
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         qicen(i,j,k) = 0.D0
         hnew (i,j,k) = 0.D0
        ENDDO
       ENDDO
      ENDDO

      DO j = jMin, jMax
       DO i = iMin, iMax
CML#ifdef ALLOW_AUTODIFF_TAMC
CML        ikey_1 = i
CML     &       + sNx*(j-1)
CML     &       + sNx*sNy*act1
CML     &       + sNx*sNy*max1*act2
CML     &       + sNx*sNy*max1*max2*act3
CML     &       + sNx*sNy*max1*max2*max3*act4
CML#endif 
CML#ifdef ALLOW_AUTODIFF_TAMC
CMLCADJ STORE frwatm(i,j) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE fzmloc(i,j) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE hice(i,j) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE hSnow1(i,j) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE icfrac(i,j) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE qic1(i,j) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE qic2(i,j) = comlev1_thsice_1, key=ikey_1
CML#endif

        IF (iceMask(i,j).GT.0.D0) THEN
         qicen(i,j,1)= qIc1(i,j)
         qicen(i,j,2)= qIc2(i,j)
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     initialize energies
         esurp(i,j) = 0.D0

c     make a local copy of evaporation
         evapLoc(i,j) = frwAtm(i,j)

C------------------------------------------------------------------------
C--   Compute growth and/or melting at the top and bottom surfaces
C------------------------------------------------------------------------

         xxMlt = aMlt*(hIce(i,j)-hMlt)
         xxFrz = aFrz*(hIce(i,j)-hFrz)
c--
         IF ( powerLawExp2 .GE. 1 ) THEN
          tmpMlt = xxMlt
          tmpFrz = xxFrz
          DO k=1,powerLawExp2
           tmpMlt = tmpMlt*tmpMlt
           tmpFrz = tmpFrz*tmpFrz
          ENDDO
          xxMlt = xxMlt*tmpMlt
          xxFrz = xxFrz*tmpFrz
          xxMlt = fracEnMelt -xxMlt
          xxFrz = fracEnFreez-xxFrz
         ENDIF
         enFrcMlt(i,j) = MAX( 0.D0, MIN( xxMlt, 1.D0 ) )
         enFrcFrz(i,j) = MAX( 0.D0, MIN( xxFrz, 1.D0 ) )

         IF (fzMlOc(i,j).GE. 0.D0) THEN
C     !-----------------------------------------------------------------
C     ! freezing conditions
C     !-----------------------------------------------------------------
          Fbot(i,j) = fzMlOc(i,j)
          IF ( icFrac(i,j).LT.iceMaskMax ) THEN
           Fbot(i,j) = enFrcFrz(i,j)*fzMlOc(i,j)
          ENDIF
         ELSE
C     !-----------------------------------------------------------------
C     ! melting conditions
C     !-----------------------------------------------------------------
C     for no currents:
          ustar = 5.D-2
C frictional velocity between ice and water
          IF (v2oc(i,j) .NE. 0.)
     &     ustar = SQRT(0.00536D0*v2oc(i,j))
          ustar=max(5.D-3,ustar)
          cpchr =cpWater*rhosw*bMeltCoef
          Fbot(i,j) = cpchr*(tFrz(i,j)-tOce(i,j))*ustar
C     fzMlOc < Fbot < 0
          Fbot(i,j) = max(Fbot(i,j),fzMlOc(i,j))
          Fbot(i,j) = min(Fbot(i,j),0.D0)
         ENDIF

C  mass of fresh water and salt initially present in ice
         mwater0(i,j) = rhos*hSnow1(i,j) + rhoi*hIce(i,j)
         msalt0 (i,j) = rhoi*hIce(i,j)*saltIce

         IF (dBug(i,j,bi,bj) ) WRITE(6,1020)
     &        'ThSI_CALC_TH: evpAtm, fzMlOc, Fbot =',
     &        frwAtm(i,j),fzMlOc(i,j),Fbot(i,j)
C     endif iceMask > 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0) THEN

C     Compute energy available for melting/growth.

         IF ( fracEnMelt.EQ.0.D0 ) THEN
          frace = 0.D0
         ELSE
          frace = (icFrac(i,j) - lowIcFrac1)/(lowIcFrac2-iceMaskMin)
          frace = MIN( enFrcMlt(i,j), MAX( 0.D0, frace ) )
         ENDIF

c     IF (tSrf(i,j) .EQ. 0.D0 .AND. sHeat(i,j).GT.0.D0) THEN
         IF ( sHeat(i,j).GT.0.D0 ) THEN
          etop(i,j) = (1.D0-frace)*sHeat(i,j) * dt
          etope(i,j) = frace*sHeat(i,j) * dt
         ELSE
          etop(i,j) =  0.D0
          etope(i,j) = 0.D0
C jmc: found few cases where tSrf=0 & sHeat < 0 : add this line to conserv energy:
          esurp(i,j) = sHeat(i,j) * dt
         ENDIF
C--   flux at the base of sea-ice:
C     conduction H.flx= flxCnB (+ =down); oceanic turbulent H.flx= Fbot (+ =down).
C-    ==> energy available(+ => melt)= (flxCnB-Fbot)*dt
c     IF (fzMlOc(i,j).LT.0.D0) THEN
c         ebot(i,j) = (1.D0-frace)*(flxCnB-Fbot(i,j)) * dt
c         ebote(i,j) = frace*(flxCnB-Fbot(i,j)) * dt
c     ELSE
c         ebot(i,j) = (flxCnB-Fbot(i,j)) * dt
c         ebote(i,j) = 0.D0
c     ENDIF
C- original formulation(above): Loose energy when flxCnB < Fbot < 0
         ebot(i,j) = (flxCnB(i,j)-Fbot(i,j)) * dt
         IF (ebot(i,j).GT.0.D0) THEN
          ebote(i,j) = frace*ebot(i,j)
          ebot(i,j)  = ebot(i,j)-ebote(i,j)
         ELSE
          ebote(i,j) = 0.D0
         ENDIF
         IF (dBug(i,j,bi,bj) ) WRITE(6,1020)
     &        'ThSI_CALC_TH: etop,etope,ebot,ebote=',
     &        etop(i,j),etope(i,j),ebot(i,j),ebote(i,j)
C     endif iceMask > 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO

C     Initialize layer thicknesses. Divide total thickness equally between
C     layers
      DO k = 1, nlyr
       DO j = jMin, jMax
        DO i = iMin, iMax
         hnew(i,j,k) = hIce(i,j) * rec_nlyr
        ENDDO
       ENDDO
      ENDDO

      DO j = jMin, jMax
       DO i = iMin, iMax
CML#ifdef ALLOW_AUTODIFF_TAMC
CMLCADJ STORE etop(i,j) = comlev1_thsice_1, key=ikey_1
CML#endif
        IF (iceMask(i,j) .GT. 0.D0 .AND.
     &         etop(i,j) .GT. 0.D0 .AND.
     &        hSnow1(i,j) .GT. 0.D0) THEN

C     Make sure internal ice temperatures do not exceed Tmlt.
C     If they do, then eliminate the layer.  (Dont think this will happen
C     for reasonable values of i0.)
C     Top melt: snow, then ice.
         rq =  rhos * qsnow
         rqh = rq * hSnow1(i,j)
         IF (etop(i,j) .LT. rqh) THEN
          hSnow1(i,j) = hSnow1(i,j) - etop(i,j)/rq
          etop(i,j) = 0.D0
         ELSE
          etop(i,j) = etop(i,j) - rqh
          hSnow1(i,j) = 0.D0
         ENDIF
C     endif iceMask > 0, etc.
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO
C     two layers of ice
      DO k = 1, nlyr
       DO j = jMin, jMax
        DO i = iMin, iMax
         IF (iceMask(i,j).GT.0.D0) THEN
CML#ifdef ALLOW_AUTODIFF_TAMC
CML          ikey_2 = k
CML     &         + nlyr*(i-1)
CML     &         + nlyr*sNx*(j-1)
CML     &         + nlyr*sNx*sNy*act1
CML     &         + nlyr*sNx*sNy*max1*act2
CML     &         + nlyr*sNx*sNy*max1*max2*act3
CML     &         + nlyr*sNx*sNy*max1*max2*max3*act4
CML#endif
CML#ifdef ALLOW_AUTODIFF_TAMC
CMLCADJ STORE etop(i,j) = comlev1_thsice_2, key=ikey_2
CMLCADJ STORE hnew(i,j,k) = comlev1_thsice_2, key=ikey_2
CML#endif
          IF (etop(i,j) .GT. 0.D0) THEN
           rq =  rhoi * qicen(i,j,k)
           rqh = rq * hnew(i,j,k)
           IF (etop(i,j) .LT. rqh) THEN
            hnew(i,j,k) = hnew(i,j,k) - etop(i,j) / rq
            etop(i,j) = 0.D0
           ELSE
            etop(i,j) = etop(i,j) - rqh
            hnew(i,j,k) = 0.D0
           ENDIF
          ELSE
           etop(i,j)=0.D0
          ENDIF
C If ice is gone and melting energy remains
c     IF (etop(i,j) .GT. 0.D0) THEN
c        WRITE (6,*)  'QQ All ice melts from top  ', i,j
c        hIce(i,j)=0.D0
c        go to 200
c     ENDIF

C     endif iceMask > 0
         ENDIF
C     end i/j-loops
        ENDDO
       ENDDO
C     end k-loop
      ENDDO

C Bottom growth:
      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0 .AND. ebot(i,j) .LT. 0.D0) THEN
C Compute enthalpy of new ice growing at bottom surface.
         qbot =  -cpIce *tFrz(i,j) + Lfresh
         dhi = -ebot(i,j) / (qbot * rhoi)
         ebot(i,j) = 0.D0
cph         k = nlyr
CML#ifdef ALLOW_AUTODIFF_TAMC
CMLCADJ STORE hnew(i,j,:) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE qicen(i,j,:) = comlev1_thsice_1, key=ikey_1
CML#endif
         qicen(i,j,nlyr) =
     &        (hnew(i,j,nlyr)*qicen(i,j,nlyr)+dhi*qbot) /
     &        (hnew(i,j,nlyr)+dhi)
         hnew(i,j,nlyr) = hnew(i,j,nlyr) + dhi
         frzSeaWat(i,j) = rhoi*dhi/dt

C     endif iceMask > 0 and ebot < 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO


C Bottom melt:
      DO k = nlyr, 1, -1
CML#ifdef ALLOW_AUTODIFF_TAMC
CML           ikey_2 = (nlyr-k+1)
CML     &         + nlyr*(i-1)
CML     &         + nlyr*sNx*(j-1)
CML     &         + nlyr*sNx*sNy*act1
CML     &         + nlyr*sNx*sNy*max1*act2
CML     &         + nlyr*sNx*sNy*max1*max2*act3
CML     &         + nlyr*sNx*sNy*max1*max2*max3*act4
CML#endif
CML#ifdef ALLOW_AUTODIFF_TAMC
CMLCADJ STORE ebot(i,j) = comlev1_thsice_2, key=ikey_2
CMLCADJ STORE hnew(i,j,k) = comlev1_thsice_2, key=ikey_2
CMLCADJ STORE qicen(i,j,k) = comlev1_thsice_2, key=ikey_2
CML#endif
       DO j = jMin, jMax
        DO i = iMin, iMax
         IF (iceMask(i,j) .GT. 0.D0 .AND.
     &        ebot(i,j)   .GT. 0.D0 .AND.
     &        hnew(i,j,k) .GT. 0.D0) THEN
          rq =  rhoi * qicen(i,j,k)
          rqh = rq * hnew(i,j,k)
          IF (ebot(i,j) .LT. rqh) THEN
           hnew(i,j,k) = hnew(i,j,k) - ebot(i,j) / rq
           ebot(i,j) = 0.D0
          ELSE
           ebot(i,j) = ebot(i,j) - rqh
           hnew(i,j,k) = 0.D0
          ENDIF
C     endif iceMask > 0 etc.
         ENDIF
C     end i/j-loops
        ENDDO
       ENDDO
C     end k-loop
      ENDDO
C     If ice melts completely and snow is left, remove the snow with
C     energy from the mixed layer
      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j) .GT. 0.D0 .AND.
     &       ebot(i,j)   .GT. 0.D0 .AND.
     &       hSnow1(i,j)  .GT. 0.D0) THEN
         rq =  rhos * qsnow
         rqh = rq * hSnow1(i,j)
         IF (ebot(i,j) .LT. rqh) THEN
          hSnow1(i,j) = hSnow1(i,j) - ebot(i,j) / rq
          ebot(i,j) = 0.D0
         ELSE
          ebot(i,j) = ebot(i,j) - rqh
          hSnow1(i,j) = 0.D0
         ENDIF
c        IF (ebot(i,j) .GT. 0.D0) THEN
c           IF (dBug) WRITE(6,*) 'All ice (& snow) melts from bottom'
c           hIce(i,j)=0.D0
c           go to 200
c        ENDIF

C     endif iceMask > 0, etc.
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO
      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0) THEN
C Compute new total ice thickness.
         hIce(i,j) = hnew(i,j,1) + hnew(i,j,2)
         IF (dBug(i,j,bi,bj) ) WRITE(6,1020)
     &        'ThSI_CALC_TH:   etop, ebot, hIce, hSnow1 =',
     &        etop(i,j), ebot(i,j), hIce(i,j), hSnow1(i,j)

C If hIce < hIceMin, melt the ice.
         IF ( hIce(i,j).LT.hIceMin
     &        .AND. (hIce(i,j)+hSnow1(i,j)).GT.0.D0 ) THEN
          esurp(i,j) = esurp(i,j) - rhos*qsnow*hSnow1(i,j)
     &         - rhoi*qicen(i,j,1)*hnew(i,j,1)
     &         - rhoi*qicen(i,j,2)*hnew(i,j,2)
          hIce(i,j)   = 0.D0
          hSnow1(i,j)  = 0.D0
          tSrf(i,j)   = 0.D0
          icFrac(i,j) = 0.D0
          qicen(i,j,1) = 0.D0
          qicen(i,j,2) = 0.D0
          IF (dBug(i,j,bi,bj) ) WRITE(6,1020)
     &         'ThSI_CALC_TH: -1 : esurp=',esurp(i,j)
         ENDIF

C     endif iceMask > 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO

      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0) THEN

C--   do a mass-budget of sea-ice to compute "fresh" = the fresh-water flux
C     that is returned to the ocean ; needs to be done before snow/evap
         frw2oc(i,j) = (mwater0(i,j)
     &        - (rhos*hSnow1(i,j)+rhoi*hIce(i,j)))/dt

         IF ( hIce(i,j) .LE. 0.D0 ) THEN
C-    return  snow to the ocean (account for Latent heat of freezing)
          frw2oc(i,j) = frw2oc(i,j) + snowP(i,j)
          flx2oc(i,j) = flx2oc(i,j) - snowP(i,j)*Lfresh
         ENDIF

C     endif iceMask > 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO
C-    else: hIce > 0
      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0) THEN

         IF ( hIce(i,j) .GT. 0.D0 ) THEN
C Let it snow
          hSnow1(i,j) = hSnow1(i,j) + dt*snowP(i,j)/rhos
C If ice evap is used to sublimate surface snow/ice or
C if no ice pass on to ocean
CML#ifdef ALLOW_AUTODIFF_TAMC
CMLCADJ STORE evapLoc(i,j) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE hSnow1(i,j) = comlev1_thsice_1, key=ikey_1
CML#endif
          IF (hSnow1(i,j).GT.0.D0) THEN
           IF (evapLoc(i,j)/rhos *dt.GT.hSnow1(i,j)) THEN
            evapLoc(i,j)=evapLoc(i,j)-hSnow1(i,j)*rhos/dt
            hSnow1(i,j)=0.D0
           ELSE
            hSnow1(i,j) = hSnow1(i,j) - evapLoc(i,j)/rhos *dt
            evapLoc(i,j)=0.D0
           ENDIF
          ENDIF
C     endif hice > 0
         ENDIF
C     endif iceMask > 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO


C-    else: hIce > 0
      DO k = 1, nlyr
       DO j = jMin, jMax
        DO i = iMin, iMax
         IF (iceMask(i,j).GT.0.D0 ) THEN

CML#ifdef ALLOW_AUTODIFF_TAMC
CMLCADJ STORE evapLoc(i,j) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE hIce(i,j) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE hnew(i,j,:) = comlev1_thsice_1, key=ikey_1
CMLCADJ STORE qicen(i,j,:) = comlev1_thsice_1, key=ikey_1
CML#endif
           IF (hIce(i,j).GT.0.D0.AND.evapLoc(i,j).GT.0.D0) THEN
CML#ifdef ALLOW_AUTODIFF_TAMC
CML            ikey_2 = k
CML     &         + nlyr*(i-1)
CML     &         + nlyr*sNx*(j-1)
CML     &         + nlyr*sNx*sNy*act1
CML     &         + nlyr*sNx*sNy*max1*act2
CML     &         + nlyr*sNx*sNy*max1*max2*act3
CML     &         + nlyr*sNx*sNy*max1*max2*max3*act4
CML#endif
CMLC--
CML#ifdef ALLOW_AUTODIFF_TAMC
CMLCADJ STORE evapLoc(i,j) = comlev1_thsice_2, key=ikey_2
CMLCADJ STORE hnew(i,j,k) = comlev1_thsice_2, key=ikey_2
CMLCADJ STORE qicen(i,j,k) = comlev1_thsice_2, key=ikey_2
CML#endif
C            IF (evapLoc(i,j) .GT. 0.D0) THEN
C-- original scheme, does not care about ice temp.
C-  this can produce small error (< 1.W/m2) in the Energy budget
c              IF (evapLoc(i,j)/rhoi *dt.GT.hnew(i,j,k)) THEN
c                evapLoc(i,j)=evapLoc(i,j)-hnew(i,j,k)*rhoi/dt
c                hnew(i,j,k)=0.D0
c              ELSE
c                hnew(i,j,k) = hnew(i,j,k) - evapLoc(i,j)/rhoi *dt
c                evapLoc(i,j)=0.D0
c              ENDIF
C-- modified scheme. taking into account Ice enthalpy
             dhi = evapLoc(i,j)/rhoi*dt
             IF (dhi.GE.hnew(i,j,k)) THEN
              evapLoc(i,j)=evapLoc(i,j)-hnew(i,j,k)*rhoi/dt
              esurp(i,j) = esurp(i,j)
     &             - hnew(i,j,k)*rhoi*(qicen(i,j,k)-Lfresh)
              hnew(i,j,k)=0.D0
             ELSE
CML#ifdef ALLOW_AUTODIFF_TAMC
CMLCADJ STORE hnew(i,j,k) = comlev1_thsice_2, key=ikey_2
CML#endif
              hq = hnew(i,j,k)*qicen(i,j,k)-dhi*Lfresh
              hnew(i,j,k) = hnew(i,j,k) - dhi
              qicen(i,j,k)=hq/hnew(i,j,k)
              evapLoc(i,j)=0.D0
             ENDIF
C-------
c     IF (evapLoc(i,j) .GT. 0.D0) THEN
c           WRITE (6,*)  'BB All ice sublimates', i,j
c           hIce(i,j)=0.D0
c           go to 200
c     ENDIF
C     endif hice > 0 and evaploc > 0
          ENDIF
C     endif iceMask > 0
         ENDIF
C     end i/j-loops
        ENDDO
       ENDDO
C     end k-loop
      ENDDO


C     still else: hice > 0
      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0) THEN
         IF (hIce(i,j) .GT. 0.D0) THEN
C Compute new total ice thickness.
          hIce(i,j) = hnew(i,j,1) + hnew(i,j,2)
C If hIce < hIceMin, melt the ice.
          IF ( hIce(i,j).GT.0.D0 .AND. hIce(i,j).LT.hIceMin ) THEN
           frw2oc(i,j) = frw2oc(i,j)
     &          + (rhos*hSnow1(i,j) + rhoi*hIce(i,j))/dt
           esurp(i,j) = esurp(i,j) - rhos*qsnow*hSnow1(i,j)
     &          - rhoi*qicen(i,j,1)*hnew(i,j,1)
     &          - rhoi*qicen(i,j,2)*hnew(i,j,2)
           hIce(i,j)   = 0.D0
           hSnow1(i,j)  = 0.D0
           tSrf(i,j)   = 0.D0
           icFrac(i,j) = 0.D0
           qicen(i,j,1) = 0.D0
           qicen(i,j,2) = 0.D0
           IF (dBug(i,j,bi,bj) ) WRITE(6,1020)
     &          'ThSI_CALC_TH: -2 : esurp,frw2oc=',
     &          esurp(i,j), frw2oc(i,j)
          ENDIF

C--   else hIce > 0: end
         ENDIF

C     endif iceMask > 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO


      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0) THEN

         IF ( hIce(i,j) .GT. 0.D0 ) THEN

C If there is enough snow to lower the ice/snow interface below
C freeboard, convert enough snow to ice to bring the interface back
C to sea-level OR if snow height is larger than hsMax, snow is
C converted to ice to bring hSnow1 down to hsMax. Largest change is
C applied and enthalpy of top ice layer adjusted accordingly.


          IF ( hSnow1(i,j) .GT. hIce(i,j)*floodFac
     &         .OR. hSnow1(i,j) .GT. hsMax ) THEN
cBB               WRITE (6,*)  'Freeboard adjusts'
c          dhi = (hSnow1(i,j) * rhos - hIce(i,j) * rhoiw) / rhosw
c          dhs = dhi * rhoi / rhos
           dhs = (hSnow1(i,j) - hIce(i,j)*floodFac) * rhoi / rhosw
           dhs = MAX( hSnow1(i,j) - hsMax, dhs )
           dhi = dhs * rhos / rhoi
           rqh = rhoi*qicen(i,j,1)*hnew(i,j,1) + rhos*qsnow*dhs
           hnew(i,j,1)    = hnew(i,j,1) + dhi
           qicen(i,j,1)   = rqh / (rhoi*hnew(i,j,1))
           hIce(i,j)  = hIce(i,j) + dhi
           hSnow1(i,j) = hSnow1(i,j) - dhs
          ENDIF

C limit ice height
C- NOTE: this part does not conserve Energy ;
C        but surplus of fresh water and salt are taken into account.
          IF (hIce(i,j).GT.hiMax) THEN
cBB      print*,'BBerr, hIce>hiMax',i,j,hIce(i,j)
           chi=hIce(i,j)-hiMax
           hnew(i,j,1)=hnew(i,j,1)-chi/2.D0
           hnew(i,j,2)=hnew(i,j,2)-chi/2.D0
           frw2oc(i,j) = frw2oc(i,j) + chi*rhoi/dt
          ENDIF
c       IF (hSnow1(i,j).GT.hsMax) THEN
cc        print*,'BBerr, hSnow1>hsMax',i,j,hSnow1(i,j)
c         chs=hSnow1(i,j)-hsMax
c         hSnow1(i,j)=hsMax
c         frw2oc(i,j) = frw2oc(i,j) + chs*rhos/dt
c       ENDIF

C Compute new total ice thickness.
          hIce(i,j) = hnew(i,j,1) + hnew(i,j,2)

          IF (dBug(i,j,bi,bj) ) WRITE(6,1020)
     &         'ThSI_CALC_TH: b-Winton: hnew, qice =',
     &         hnew(i,j,1), hnew(i,j,2),
     &         qicen(i,j,1), qicen(i,j,2)

          hlyr = hIce(i,j) * rec_nlyr
CML          CALL THSICE_RESHAPE_LAYERS(
CML     U         qicen(i,j,:),
CML     I         hlyr, hnew(i,j,:), myThid )
C     inlined version of S/R THSICE_RESHAPE_LAYERS
C     | Repartition into equal-thickness layers, conserving energy.
C     *==========================================================*
C     | This is the 2-layer version (formerly "NEW_LAYERS_WINTON")
C     |  from M. Winton 1999, JAOT, sea-ice model.
          if (hnew(i,j,1).gt.hnew(i,j,2)) then
C--   Layer 1 gives ice to layer 2
           f1 = (hnew(i,j,1)-hlyr)/hlyr
           q2tmp = f1*qicen(i,j,1) + (1.D0-f1)*qicen(i,j,2)
           if (q2tmp.gt.Lfresh) then
            qicen(i,j,2) = q2tmp
           else
C-    Keep q2 fixed to avoid q2<Lfresh and T2>0
            qh2 = hlyr*qicen(i,j,2)
            qhtot = hnew(i,j,1)*qicen(i,j,1) + hnew(i,j,2)*qicen(i,j,2)
            qh1 = qhtot - qh2
            qicen(i,j,1) = qh1/hlyr
           endif
          else
C-    Layer 2 gives ice to layer 1
           f1 = hnew(i,j,1)/hlyr
           qicen(i,j,1) = f1*qicen(i,j,1) + (1.D0-f1)*qicen(i,j,2)
          endif
C     end of inlined S/R THSICE_RESHAPE_LAYERS

          IF (dBug(i,j,bi,bj) ) WRITE(6,1020)
     &         'ThSI_CALC_TH: icFrac,hIce, qtot, hSnow1 =',
     &         icFrac(i,j),hIce(i,j), (qicen(i,j,1)+qicen(i,j,2))*0.5,
     &         hSnow1(i,j)

C-    if hIce > 0 : end
         ENDIF
c200     CONTINUE

C     endif iceMask > 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO

      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0) THEN
C-  Compute surplus energy left over from melting.

         IF (hIce(i,j).LE.0.D0) icFrac(i,j)=0.D0

C.. heat fluxes left over for ocean
         flx2oc(i,j) = flx2oc(i,j)
     &        + (Fbot(i,j)+(esurp(i,j)+etop(i,j)+ebot(i,j))/dt)
         IF (dBug(i,j,bi,bj) ) WRITE(6,1020)
     &        'ThSI_CALC_TH: [esurp,etop+ebot]/dt =',
     &        esurp(i,j)/dt,etop(i,j)/dt,ebot(i,j)/dt

C--   Evaporation left to the ocean :
         frw2oc(i,j) = frw2oc(i,j) - evapLoc(i,j)
C--   Correct Atmos. fluxes for this different latent heat:
C     evap was computed over freezing surf.(tSrf<0), latent heat = Lvap+Lfresh
C     but should be Lvap only for the fraction "evap" that is left to the ocean.
         flx2oc(i,j) = flx2oc(i,j) + evapLoc(i,j)*Lfresh

C fresh and salt fluxes
c     frw2oc(i,j) = (mwater0(i,j) - (rhos*(hSnow1(i,j))
c    &              + rhoi*(hIce(i,j))))/dt-evapLoc(i,j)
c     fsalt = (msalt0(i,j) - rhoi*hIce(i,j)*saltIce)/35.D0/dt  ! for same units as frw2oc
C     note (jmc): frw2oc is computed from a sea-ice mass budget that already
C     contains, at this point, snow & evaporation (of snow & ice)
C     but are not meant to be part of ice/ocean fresh-water flux.
C     fix: a) like below or b) by making the budget before snow/evap is added
c     frw2oc(i,j) = (mwater0(i,j) - (rhos*(hSnow1(i,j)) + rhoi*(hIce(i,j))))/dt
c    &      + snow(i,j,bi,bj)*rhos - frwAtm
         fsalt(i,j) = (msalt0(i,j) - rhoi*hIce(i,j)*saltIce)/dt

         IF (dBug(i,j,bi,bj) ) THEN
          WRITE(6,1020)'ThSI_CALC_TH:dH2O,Ev[kg],frw2oc,fsalt',
     &         (mwater0(i,j)-(rhos*hSnow1(i,j)+rhoi*hIce(i,j)))/dt,
     &         evapLoc(i,j),frw2oc(i,j),fsalt(i,j)
          WRITE(6,1020)'ThSI_CALC_TH: flx2oc,Fbot,extend/dt =',
     &         flx2oc(I,J),Fbot(i,j),(etope(i,j)+ebote(i,j))/dt
         ENDIF

C--   add remaining liquid Precip (rain+RunOff) directly to ocean:
         frw2oc(i,j) = frw2oc(i,j) + (prcAtm(i,j)-snowP(i,j))

C     endif iceMask > 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO


      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0) THEN
C--   note: at this point, icFrac has not been changed (unless reset to zero)
C     and it can only be reduced by lateral melting in the following part:

C     calculate extent changes
         extend=etope(i,j)+ebote(i,j)
         IF (icFrac(i,j).GT.0.D0.AND.extend.GT.0.D0) THEN
          rq =  rhoi * 0.5D0*(qicen(i,j,1)+qicen(i,j,2))
          rs =  rhos * qsnow
          rqh = rq * hIce(i,j) + rs * hSnow1(i,j)
          freshe=(rhos*hSnow1(i,j)+rhoi*hIce(i,j))/dt
          salte=(rhoi*hIce(i,j)*saltIce)/dt
          IF ( extend.LT.rqh ) THEN
           icFrac(i,j)=(1.D0-extend/rqh)*icFrac(i,j)
          ENDIF
          IF ( extend.LT.rqh .AND. icFrac(i,j).GE.iceMaskMin ) THEN
           frw2oc(i,j)=frw2oc(i,j)+extend/rqh*freshe
           fsalt(i,j)=fsalt(i,j)+extend/rqh*salte
          ELSE
           icFrac(i,j)=0.D0
           hIce(i,j)  =0.D0
           hSnow1(i,j) =0.D0
           flx2oc(i,j)=flx2oc(i,j)+(extend-rqh)/dt
           frw2oc(i,j)=frw2oc(i,j)+freshe
           fsalt(i,j)=fsalt(i,j)+salte
          ENDIF
         ELSEIF (extend.GT.0.D0) THEN
          flx2oc(i,j)=flx2oc(i,j)+extend/dt
         ENDIF
C     endif iceMask > 0
        ENDIF
C     end i/j-loops
       ENDDO
      ENDDO
      DO j = jMin, jMax
       DO i = iMin, iMax
        IF (iceMask(i,j).GT.0.D0) THEN
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Update output variables :

C--   Diagnostic of Atmos. fresh water flux (E-P) over sea-ice :
C     substract precip from Evap (<- stored in frwAtm array)
         frwAtm(i,j) = frwAtm(i,j) - prcAtm(i,j)

C--   update Mixed-Layer Freezing potential heat flux by substracting the
C     part which has already been accounted for (Fbot):
         fzMlOc(i,j) = fzMlOc(i,j) - Fbot(i,j)*iceMask(i,j)

C-- Update Sea-Ice state output:
         qIc1(i,j)   = qicen(i,j,1)
         qIc2(i,j)   = qicen(i,j,2)
         IF (dBug(i,j,bi,bj) ) WRITE(6,1020)
     &        'ThSI_CALC_TH: icFrac,flx2oc,fsalt,frw2oc=',
     &        icFrac(i,j), flx2oc(i,j), fsalt(i,j), frw2oc(i,j)
C     endif iceMask > 0
        ENDIF
       ENDDO
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      RETURN
      END
