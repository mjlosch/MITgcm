



















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
C     !ROUTINE: THSICE_RESHAPE_LAYERS
C     !INTERFACE:
      SUBROUTINE THSICE_RESHAPE_LAYERS(
     U                                 qicen,
     I                                 hlyr, hnew, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R THSICE_RESHAPE_LAYERS
C     | Repartition into equal-thickness layers, conserving energy.
C     *==========================================================*
C     | This is the 2-layer version (formerly "NEW_LAYERS_WINTON") 
C     |  from M. Winton 1999, JAOT, sea-ice model.
C     *==========================================================*
C     \ev

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
C     qicen  :: ice enthalpy (J/kg)
C     hnew   :: new ice layer thickness (m)
C     hlyr   :: individual ice layer thickness (m)
C     myThid :: Number of this instance
      Real*8  qicen(*)
      Real*8  hnew(*)
      Real*8  hlyr
      INTEGER myThid
CEOP

C     == Local Variables ==
C      f1         :: Fraction of upper layer ice in new layer
C      qh1, qh2   :: qice*h for layers 1 and 2
C      qhtot      :: qh1 + qh2
C      q2tmp      :: Temporary value of qice for layer 2
      Real*8  f1
      Real*8  qh1, qh2
      Real*8  qhtot
      Real*8  q2tmp

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      if (hnew(1).gt.hnew(2)) then
C-    Layer 1 gives ice to layer 2
         f1 = (hnew(1)-hlyr)/hlyr
         q2tmp = f1*qicen(1) + (1.D0-f1)*qicen(2)
         if (q2tmp.gt.Lfresh) then
            qicen(2) = q2tmp
         else
C-    Keep q2 fixed to avoid q2<Lfresh and T2>0
            qh2 = hlyr*qicen(2)
            qhtot = hnew(1)*qicen(1) + hnew(2)*qicen(2)
            qh1 = qhtot - qh2
            qicen(1) = qh1/hlyr
         endif
      else
C-    Layer 2 gives ice to layer 1
         f1 = hnew(1)/hlyr
         qicen(1) = f1*qicen(1) + (1.D0-f1)*qicen(2)
      endif


      RETURN
      END
