



















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
C     !ROUTINE: THSICE_SOLVE4TEMP
C     !INTERFACE:
      SUBROUTINE THSICE_SOLVE4TEMP(
     I                  bi, bj,
     I                  iMin,iMax, jMin,jMax, dBugFlag,
     I                  useBulkForce, useEXF,
     I                  icMask, hIce, hSnow1, tFrz, flxExSW,
     U                  flxSW, tSrf1, qIc1, qIc2,
     O                  tIc1, tIc2, dTsrf1,
     O                  sHeat, flxCnB, flxAtm, evpAtm,
     I                  myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R  THSICE_SOLVE4TEMP
C     *==========================================================*
C     | Solve (implicitly) for sea-ice and surface temperature
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

C     !USES:
      IMPLICIT NONE

C     == Global variables ===
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
C     useBulkForce:: use surf. fluxes from bulk-forcing external S/R
C     useEXF      :: use surf. fluxes from exf          external S/R
C---  Input:
C     icMask      :: sea-ice fractional mask [0-1]
C     hIce        :: ice height [m]
C     hSnow1      :: snow height [m]
C     tFrz        :: sea-water freezing temperature [oC] (function of S)
C     flxExSW     :: surf. heat flux (+=down) except SW, function of surf. temp Ts:
C                    0: Flx(Ts=0) ; 1: Flx(Ts=Ts^n) ; 2: d.Flx/dTs(Ts=Ts^n)
C---  Modified (input&output):
C     flxSW       :: net Short-Wave flux (+=down) [W/m2]: input= at surface
C                 ::                 output= below sea-ice, into the ocean
C     tSrf1       :: surface (ice or snow) temperature
C     qIc1        :: ice enthalpy (J/kg), 1rst level
C     qIc2        :: ice enthalpy (J/kg), 2nd level
C---  Output
C     tIc1        :: temperature of ice layer 1 [oC]
C     tIc2        :: temperature of ice layer 2 [oC]
C     dTsrf1      :: surf. temp adjusment: Ts^n+1 - Ts^n
C     sHeat       :: surf heating flux left to melt snow or ice (= Atmos-conduction)
C     flxCnB      :: heat flux conducted through the ice to bottom surface
C     flxAtm      :: net flux of energy from the atmosphere [W/m2] (+=down)
C                    without snow precip. (energy=0 for liquid water at 0.oC)
C     evpAtm      :: evaporation to the atmosphere [kg/m2/s] (>0 if evaporate)
C---  Input:
C     myTime      :: current Time of simulation [s]
C     myIter      :: current Iteration number in simulation
C     myThid      :: my Thread Id number
      INTEGER bi,bj
      INTEGER iMin, iMax
      INTEGER jMin, jMax
      LOGICAL dBugFlag
      LOGICAL useBulkForce
      LOGICAL useEXF
      Real*8 icMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 hIce   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 hSnow1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 tFrz   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 flxSW  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 tSrf1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 qIc1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 qIc2   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 tIc1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 tIc2   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 sHeat  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 flxCnB (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 flxAtm (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 evpAtm (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 flxExSW(iMin:iMax,jMin:jMax,0:2)
      Real*8 dTsrf1 (iMin:iMax,jMin:jMax)
      Real*8 myTime
      INTEGER myIter
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     == Local Variables ==
C     useBlkFlx    :: use some bulk-formulae to compute fluxes over sea-ice
C                  :: (otherwise, fluxes are passed as argument from AIM)
C     iterate4Tsf  :: to stop to iterate when all icy grid pts Tsf did converged
C     iceFlag      :: True= do iterate for Surf.Temp ; False= do nothing
C     frsnow       :: fractional snow cover
C     fswpen       :: SW penetrating beneath surface (W m-2)
C     fswdn        :: SW absorbed at surface (W m-2)
C     fswint       :: SW absorbed in ice (W m-2)
C     fswocn       :: SW passed through ice to ocean (W m-2)
C     Tsf          :: surface (ice or snow) temperature   (local copy of tSrf1)
C     flx0exSW     :: net surface heat flux over melting snow/ice, except short-wave.
C     flxTexSW     :: net surface heat flux, except short-wave (W/m2)
C     dFlxdT       :: deriv of flxNet wrt Tsf (W m-2 deg-1)
C     evap00       :: evaporation over melting snow/ice [kg/m2/s] (>0 if evaporate)
C                     renamed to evap00 because TAF confuses symbol with EXF evap0
C     evapT        :: evaporation over snow/ice [kg/m2/s] (>0 if evaporate)
C     dEvdT        :: derivative of evap. with respect to Tsf [kg/m2/s/K]
C     flxNet       :: net surf heat flux (+=down), from Atmos. to sea-ice (W m-2)
C     netSW        :: net Short-Wave flux at surface (+=down) [W/m2]
C     fct          :: heat conducted to top surface
C     k12, k32     :: thermal conductivity terms
C     a10,b10,c10  :: coefficients in quadratic eqn for T1
C     a1, b1, c1   :: coefficients in quadratic eqn for T1
C     dt           :: timestep
      LOGICAL useBlkFlx
      LOGICAL iterate4Tsf
      INTEGER i, j, k, iterMax
      INTEGER ii, jj, icount, stdUnit
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      Real*8  frsnow
      Real*8  fswpen
      Real*8  fswdn
      Real*8  fswint
      Real*8  fswocn
      Real*8  iceFlag (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  Tsf     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  flx0exSW(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  flxTexSW(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  dFlxdT  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  evap00  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  evapT   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  dEvdT   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  flxNet
      Real*8  fct
      Real*8  k12     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  k32
      Real*8  a10     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  b10     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  c10     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8  a1, b1, c1
      Real*8  dt
      Real*8  recip_dhSnowLin
      Real*8  netSW

C-    Define grid-point location where to print debugging values
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

 1010 FORMAT(A,I3,3F11.6)
 1020 FORMAT(A,1P4E14.6)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|



      stdUnit = standardMessageUnit
      useBlkFlx = useEXF .OR. useBulkForce
      dt  = thSIce_dtTemp
      IF ( dhSnowLin.GT.0.D0 ) THEN
        recip_dhSnowLin = 1.D0 / dhSnowLin
      ELSE
        recip_dhSnowLin = 0.D0
      ENDIF

      iterate4Tsf = .FALSE.
      icount = 0

      DO j = jMin, jMax
       DO i = iMin, iMax
        IF ( icMask(i,j).GT.0.D0) THEN
          iterate4Tsf  = .TRUE.
          iceFlag(i,j) = 1.D0
          IF ( dBug(i,j,bi,bj) ) WRITE(stdUnit,'(A,2I4,2I2)')
     &         'ThSI_SOLVE4T: i,j=',i,j,bi,bj
          IF ( hIce(i,j).LT.hIceMin ) THEN
C     if hi < hIceMin, melt the ice.
C     keep the position of this problem but do the stop later
           ii = i
           jj = j
           icount = icount + 1
          ENDIF

C--   Fractional snow cover:
C     assume a linear distribution of snow thickness, with dhSnowLin slope,
C      from hs-dhSnowLin to hs+dhSnowLin if full ice & snow cover.
C     frsnow = fraction of snow over the ice-covered part of the grid cell
          IF ( hSnow1(i,j) .GT. icMask(i,j)*dhSnowLin ) THEN
           frsnow = 1.D0
          ELSE
           frsnow = hSnow1(i,j)*recip_dhSnowLin/icMask(i,j)
           IF ( frsnow.GT.0.D0 ) frsnow = SQRT(frsnow)
          ENDIF

C--   Compute SW flux absorbed at surface and penetrating to layer 1.
          fswpen = flxSW(i,j) * (1.D0 - frsnow) * i0swFrac
          fswocn = fswpen * exp(-ksolar*hIce(i,j))
          fswint = fswpen - fswocn
          fswdn  = flxSW(i,j) - fswpen

C     Initialise Atmospheric surf. heat flux with net SW flux at the surface
          flxAtm(i,j) = flxSW(i,j)
C     SW flux at sea-ice base left to the ocean
          flxSW(i,j) = fswocn
C     Initialise surface Heating with SW contribution
          sHeat(i,j) = fswdn

C--   Compute conductivity terms at layer interfaces.
          k12(i,j) = 4.D0*kIce*kSnow
     &             / (kSnow*hIce(i,j) + 4.D0*kIce*hSnow1(i,j))
          k32      = 2.D0*kIce  / hIce(i,j)

C--   Compute ice temperatures
          a1 = cpIce
          b1 = qIc1(i,j) + (cpWater-cpIce )*Tmlt1 - Lfresh
          c1 = Lfresh * Tmlt1
          tIc1(i,j) = 0.5D0 *(-b1 - SQRT(b1*b1-4.D0*a1*c1))/a1
          tIc2(i,j) = (Lfresh-qIc2(i,j)) / cpIce

          IF (tIc1(i,j).GT.0.D0 ) THEN
           WRITE(stdUnit,'(A,I12,1PE14.6)')
     &       ' BBerr: Tice(1) > 0 ; it=', myIter, qIc1(i,j)
           WRITE(stdUnit,'(A,4I5,2F11.4)')
     &      ' BBerr: i,j,bi,bj,Tice = ',i,j,bi,bj,tIc1(i,j),tIc2(i,j)
          ENDIF
          IF ( tIc2(i,j).GT.0.D0) THEN
           WRITE(stdUnit,'(A,I12,1PE14.6)')
     &       ' BBerr: Tice(2) > 0 ; it=', myIter, qIc2(i,j)
           WRITE(stdUnit,'(A,4I5,2F11.4)')
     &      ' BBerr: i,j,bi,bj,Tice = ',i,j,bi,bj,tIc1(i,j),tIc2(i,j)
          ENDIF
          IF ( dBug(i,j,bi,bj) ) WRITE(stdUnit,1010)
     &     'ThSI_SOLVE4T: k, Ts, Tice=',0,tSrf1(i,j),tIc1(i,j),tIc2(i,j)

C--   Compute coefficients used in quadratic formula.

          a10(i,j) = rhoi*cpIce *hIce(i,j)/(2.D0*dt) +
     &          k32*( 4.D0*dt*k32 + rhoi*cpIce *hIce(i,j) )
     &           / ( 6.D0*dt*k32 + rhoi*cpIce *hIce(i,j) )
          b10(i,j) = -hIce(i,j)*
     &          ( rhoi*cpIce*tIc1(i,j) + rhoi*Lfresh*Tmlt1/tIc1(i,j) )
     &           /(2.D0*dt)
     &        - k32*( 4.D0*dt*k32*tFrz(i,j)
     &               +rhoi*cpIce*hIce(i,j)*tIc2(i,j) )
     &           / ( 6.D0*dt*k32 + rhoi*cpIce *hIce(i,j) )
     &        - fswint
          c10(i,j) = rhoi*Lfresh*hIce(i,j)*Tmlt1 / (2.D0*dt)

        ELSE
          iceFlag(i,j) = 0.D0
        ENDIF
       ENDDO
      ENDDO
      IF ( icount .gt. 0 ) THEN
       WRITE(stdUnit,'(A,I5,A)')
     &      'THSICE_SOLVE4TEMP: there are ',icount,
     &      ' case(s) where hIce<hIceMin;'
       WRITE(stdUnit,'(A,I3,A1,I3,A)')
     &      'THSICE_SOLVE4TEMP: the last one was at (',ii,',',jj,
     &      ') with hIce = ', hIce(ii,jj)
       WRITE( msgBuf, '(A)')
     &      'THSICE_SOLVE4TEMP: should not enter if hIce<hIceMin'
       CALL PRINT_ERROR( msgBuf , myThid )
       STOP 'ABNORMAL END: S/R THSICE_SOLVE4TEMP'
      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|


C--   Get surface fluxes over melting surface
      IF ( useBlkFlx .AND. iterate4Tsf  ) THEN
        DO j = jMin, jMax
         DO i = iMin, iMax
           Tsf(i,j) = 0.
         ENDDO
        ENDDO
        IF ( useEXF ) THEN
           k = 1
           CALL THSICE_GET_EXF(   bi, bj, k,
     I                            iMin, iMax, jMin, jMax,
     I                            iceFlag, hSnow1, Tsf,
     O                            flx0exSW, dFlxdT, evap00, dEvdT,
     I                            myTime, myIter, myThid )
C-    could add this "ifdef" to hide THSICE_GET_BULKF from TAF
        ELSEIF ( useBulkForce ) THEN
           CALL THSICE_GET_BULKF( bi, bj,
     I                            iMin, iMax, jMin, jMax,
     I                            iceFlag, hSnow1, Tsf,
     O                            flx0exSW, dFlxdT, evap00, dEvdT,
     I                            myTime, myIter, myThid )
C--- end if: IF ( useEXF ) THEN
        ENDIF
C--- end if: IF ( useBlkFlx .AND. iterate4Tsf  ) THEN
      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|
C---  Compute new surface and internal temperatures; iterate until
C     Tsfc converges.
      DO j = jMin, jMax
        DO i = iMin, iMax
          Tsf(i,j)  = tSrf1(i,j)
          dTsrf1(i,j) = Terrmax
        ENDDO
      ENDDO
      IF ( useBlkFlx ) THEN
        iterMax = nitMaxTsf
      ELSE
        iterMax = 1
      ENDIF


C ----- begin iteration  -----
      DO k = 1,iterMax
       IF ( iterate4Tsf ) THEN
        iterate4Tsf = .FALSE.
C
        IF ( useBlkFlx ) THEN

C--   Compute top surface flux.
         IF ( useEXF ) THEN
           CALL THSICE_GET_EXF(   bi, bj, k+1,
     I                            iMin, iMax, jMin, jMax,
     I                            iceFlag, hSnow1, Tsf,
     O                            flxTexSW, dFlxdT, evapT, dEvdT,
     I                            myTime, myIter, myThid )
C-    could add this "ifdef" to hide THSICE_GET_BULKF from TAF
         ELSEIF ( useBulkForce ) THEN
           CALL THSICE_GET_BULKF( bi, bj,
     I                            iMin, iMax, jMin, jMax,
     I                            iceFlag, hSnow1, Tsf,
     O                            flxTexSW, dFlxdT, evapT, dEvdT,
     I                            myTime, myIter, myThid )
C--- end if: IF ( useEXF ) THEN
         ENDIF
        ELSE
         DO j = jMin, jMax
          DO i = iMin, iMax
           IF ( iceFlag(i,j).GT.0.D0 ) THEN
             flxTexSW(i,j) = flxExSW(i,j,1)
             dFlxdT(i,j) = flxExSW(i,j,2)
           ENDIF
          ENDDO
         ENDDO
C--- end if: IF ( useBlkFlx ) THEN
        ENDIF


C--   Compute new top layer and surface temperatures.
C     If Tsfc is computed to be > 0 C, fix Tsfc = 0 and recompute T1
C     with different coefficients.
        DO j = jMin, jMax
         DO i = iMin, iMax
          IF ( iceFlag(i,j).GT.0.D0 ) THEN
           flxNet = sHeat(i,j) + flxTexSW(i,j)
           IF ( dBug(i,j,bi,bj) ) WRITE(stdUnit,1020)
     &     'ThSI_SOLVE4T: flxNet,dFlxdT,k12,D=',
     &      flxNet, dFlxdT(i,j), k12(i,j), k12(i,j)-dFlxdT(i,j)

           a1 = a10(i,j) - k12(i,j)*dFlxdT(i,j) / (k12(i,j)-dFlxdT(i,j))
           b1 = b10(i,j) - k12(i,j)*(flxNet-dFlxdT(i,j)*Tsf(i,j))
     &                   /(k12(i,j)-dFlxdT(i,j))
           c1 = c10(i,j)
           tIc1(i,j)  = -(b1 + SQRT(b1*b1-4.D0*a1*c1))/(2.D0*a1)
           dTsrf1(i,j) = (flxNet + k12(i,j)*(tIc1(i,j)-Tsf(i,j)))
     &                   /(k12(i,j)-dFlxdT(i,j))
           Tsf(i,j) = Tsf(i,j) + dTsrf1(i,j)
C
           IF ( Tsf(i,j) .GT. 0.D0 ) THEN
            IF ( dBug(i,j,bi,bj) ) WRITE(stdUnit,1010)
     &     'ThSI_SOLVE4T: k,ts,t1,dTs=',k,Tsf(i,j),tIc1(i,j),dTsrf1(i,j)
            a1 = a10(i,j) + k12(i,j)
C      note: b1 = b10 - k12*Tf0
            b1 = b10(i,j)
            tIc1(i,j) = (-b1 - SQRT(b1*b1-4.D0*a1*c1))/(2.D0*a1)
            Tsf(i,j) = 0.D0
            IF ( useBlkFlx ) THEN
             flxTexSW(i,j) = flx0exSW(i,j)
             evapT(i,j) = evap00(i,j)
             dTsrf1(i,j) = 0.D0
            ELSE
             flxTexSW(i,j) = flxExSW(i,j,0)
             dTsrf1(i,j) = 1000.
             dFlxdT(i,j) = 0.
            ENDIF
           ENDIF

C--   Check for convergence.  If no convergence, then repeat.
           IF (ABS(dTsrf1(i,j)).GE.Terrmax) THEN
            iceFlag(i,j) = 1.D0
           ELSE
            iceFlag(i,j) = 0.D0
           ENDIF
           iterate4Tsf = iterate4Tsf .OR. (iceFlag(i,j).GT.0.D0)

C     Convergence test: Make sure Tsfc has converged, within prescribed error.
C     (Energy conservation is guaranteed within machine roundoff, even
C      if Tsfc has not converged.)
C     If no convergence, then repeat.

           IF ( dBug(i,j,bi,bj) ) WRITE(stdUnit,1010)
     &    'ThSI_SOLVE4T: k,ts,t1,dTs=', k,Tsf(i,j),tIc1(i,j),dTsrf1(i,j)
           IF ( useBlkFlx .AND. k.EQ.nitMaxTsf .AND.
     &                     (iceFlag(i,j).GT.0.D0) ) THEN
            WRITE(stdUnit,'(A,4I4,I12,F15.9)')
     &       ' BB: not converge: i,j,it,hi=',i,j,bi,bj,myIter,hIce(i,j)
            WRITE(stdUnit,*)
     &        'BB: not converge: Tsf, dTsf=', Tsf(i,j), dTsrf1(i,j)
            WRITE(stdUnit,*)
     &        'BB: not converge: flxNet,dFlxT=', flxNet, dFlxdT(i,j)
            IF ( Tsf(i,j).LT.-70.D0 ) THEN
             WRITE( msgBuf, '(A,2I4,2I3,I10,F12.3)')
     &        'THSICE_SOLVE4TEMP: Too low Tsf in', i, j, bi, bj,
     &                            myIter, Tsf(i,j)
             CALL PRINT_ERROR( msgBuf , myThid )
             STOP 'ABNORMAL END: S/R THSICE_SOLVE4TEMP'
            ENDIF
           ENDIF

          ENDIF
         ENDDO
        ENDDO
C--- end if: IF ( iterate4Tsf ) THEN
       ENDIF
C--- end loop DO k = 1,iterMax
      ENDDO
C ------ end iteration ------------

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|

      DO j = jMin, jMax
       DO i = iMin, iMax
        IF ( icMask(i,j).GT.0.D0 ) THEN

C--   Compute new bottom layer temperature.
          k32       = 2.D0*kIce  / hIce(i,j)
          tIc2(i,j) = ( 2.D0*dt*k32*(tIc1(i,j)+2.D0*tFrz(i,j))
     &                 + rhoi*cpIce*hIce(i,j)*tIc2(i,j))
     &               /(6.D0*dt*k32 + rhoi*cpIce*hIce(i,j))
          IF ( dBug(i,j,bi,bj) ) WRITE(stdUnit,1010)
     &  'ThSI_SOLVE4T: k, Ts, Tice=',k,Tsf(i,j),tIc1(i,j),tIc2(i,j)
          netSW   = flxAtm(i,j)

C--   Compute final flux values at surfaces.
          tSrf1(i,j)  = Tsf(i,j)
          fct    = k12(i,j)*(Tsf(i,j)-tIc1(i,j))
          flxCnB(i,j) = 4.D0*kIce *(tIc2(i,j)-tFrz(i,j))/hIce(i,j)
          flxNet = sHeat(i,j) + flxTexSW(i,j)
          flxNet = flxNet + dFlxdT(i,j)*dTsrf1(i,j)
          IF ( useBlkFlx ) THEN
C-    needs to update also Evap (Tsf changes) since Latent heat has been updated
            evpAtm(i,j) = evapT(i,j) + dEvdT(i,j)*dTsrf1(i,j)
          ELSE
C- WARNING: Evap & +Evap*Lfresh are missing ! (but only affects Diagnostics)
            evpAtm(i,j) = 0.
          ENDIF
C-    Update energy flux to Atmos with other than SW contributions;
C     use latent heat = Lvap (energy=0 for liq. water at 0.oC)
          flxAtm(i,j) = flxAtm(i,j) + flxTexSW(i,j)
     &                + dFlxdT(i,j)*dTsrf1(i,j) + evpAtm(i,j)*Lfresh
C-    excess of energy @ surface (used for surface melting):
          sHeat(i,j) = flxNet - fct

          IF ( dBug(i,j,bi,bj) ) WRITE(stdUnit,1020)
     &     'ThSI_SOLVE4T: flxNet,fct,Dif,flxCnB=',
     &                    flxNet,fct,flxNet-fct,flxCnB(i,j)

C--   Compute new enthalpy for each layer.
          qIc1(i,j) = -cpWater*Tmlt1 + cpIce *(Tmlt1-tIc1(i,j))
     &                + Lfresh*(1.D0-Tmlt1/tIc1(i,j))
          qIc2(i,j) = -cpIce *tIc2(i,j) + Lfresh

C--   Make sure internal ice temperatures do not exceed Tmlt.
C     (This should not happen for reasonable values of i0swFrac)
          IF (tIc1(i,j) .GE. Tmlt1) THEN
           WRITE(stdUnit,'(A,2I4,2I3,1P2E14.6)')
     &     ' BBerr - Bug: IceT(1) > Tmlt',i,j,bi,bj,tIc1(i,j),Tmlt1
          ENDIF
          IF (tIc2(i,j) .GE. 0.D0) THEN
           WRITE(stdUnit,'(A,2I4,2I3,1P2E14.6)')
     &     ' BBerr - Bug: IceT(2) > 0',i,j,bi,bj,tIc2(i,j)
          ENDIF

          IF ( dBug(i,j,bi,bj) ) THEN
           WRITE(stdUnit,1020) 'ThSI_SOLV_4T: Tsf, Tice(1,2), dTsurf=',
     &           Tsf(i,j), tIc1(i,j), tIc2(i,j), dTsrf1(i,j)
           WRITE(stdUnit,1020) 'ThSI_SOLV_4T: sHeat, flxCndBt, Qice =',
     &           sHeat(i,j), flxCnB(i,j), qIc1(i,j), qIc2(i,j)
           WRITE(stdUnit,1020) 'ThSI_SOLV_4T: flxA, evpA, fxSW_bf,af=',
     &           flxAtm(i,j), evpAtm(i,j), netSW, flxSW(i,j)
          ENDIF

        ELSE
C--     ice-free grid point:
c         tIc1  (i,j) = 0.D0
c         tIc2  (i,j) = 0.D0
          dTsrf1(i,j) = 0.D0
c         sHeat (i,j) = 0.D0
c         flxCnB(i,j) = 0.D0
c         flxAtm(i,j) = 0.D0
c         evpAtm(i,j) = 0.D0

        ENDIF
       ENDDO
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|

      RETURN
      END
