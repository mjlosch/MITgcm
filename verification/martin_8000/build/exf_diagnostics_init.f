












CBOP
C !ROUTINE: EXF_OPTIONS.h
C !INTERFACE:
C #include "EXF_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for EXternal Forcing (EXF) package:
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP








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



C-- Package-specific Options & Macros go here

C   --------------------
C   pkg/exf CPP options:
C   (see also table below on how to combine options)

C   > ( EXF_VERBOSE ) < replaced with run-time integer parameter "exf_debugLev"
C
C   >>> ALLOW_ATM_WIND <<<
C       If defined, set default value of run-time param. "useAtmWind" to True.
C       If useAtmWind=True, read-in and use wind vector (uwind/vwind)
C       to compute surface wind stress.
C
C   >>> ALLOW_ATM_TEMP <<<
C       This is the main EXF option controlling air-sea buoyancy fluxes:
C      If undefined, net heat flux (Qnet) and net fresh water flux
C       (EmP or EmPmR) are set according to hfluxfile & sfluxfile setting.
C      If defined, net heat flux and net fresh water flux are computed
C       from sum of various components (radiative SW,LW + turbulent heat
C       fluxes SH,LH ; Evap, Precip and optionally RunOff) thus ignoring
C       hfluxfile & sfluxfile.
C      In addition, it allows to read-in from files atmospheric temperature
C       and specific humidity, net radiative fluxes, and precip.
C       Also enable to read-in Evap (if EXF_READ_EVAP is defined) or
C       turbulent heat fluxes (if ALLOW_READ_TURBFLUXES is defined).
C
C   >>> ALLOW_DOWNWARD_RADIATION <<<
C       If defined, downward long-wave and short-wave radiation
C       can be read-in form files to compute net lwflux and swflux.
C
C   >>> ALLOW_ZENITHANGLE <<<
C       If defined, ocean albedo varies with the zenith angle, and
C       incoming fluxes at the top of the atmosphere are computed
C
C   >>> ALLOW_BULKFORMULAE <<<
C       Allows the use of bulk formulae in order to estimate
C       turbulent fluxes (Sensible,Latent,Evap) at the ocean surface.
C
C   >>> EXF_CALC_ATMRHO
C       Calculate the local air density as function of temp, humidity
C       and pressure
C
C   >>> EXF_READ_EVAP <<<
C       If defined, evaporation field is read-in from file;
C     Note: if ALLOW_BULKFORMULAE is defined, evap that is computed from
C       atmospheric state will be replaced by read-in evap but computed
C       latent heat flux will be kept.
C
C   >>> ALLOW_READ_TURBFLUXES <<<
C       If defined, turbulent heat fluxes (sensible and latent) can be read-in
C       from files (but overwritten if ALLOW_BULKFORMULAE is defined).
C
C   >>> ALLOW_RUNOFF <<<
C       If defined, river and glacier runoff can be read-in from files.
C
C   >>> ALLOW_SALTFLX <<<
C       If defined, upward salt flux can be read-in from files.
C
C   >>> ALLOW_RUNOFTEMP <<<
C       If defined, river and glacier runoff temperature
C       can be read-in from files.
C
C   >>>  <<<
C       If defined, atmospheric pressure can be read-in from files.
C   WARNING: this flag is set (define/undef) in CPP_OPTIONS.h
C            and cannot be changed here (in EXF_OPTIONS.h)
C
C   >>> EXF_ALLOW_TIDES <<<
C       If defined, 2-D tidal geopotential can be read-in from files
C
C   >>> EXF_SEAICE_FRACTION <<<
C       If defined, seaice fraction can be read-in from files (areaMaskFile)
C
C   >>> ALLOW_CLIMSST_RELAXATION <<<
C       Allow the relaxation to a monthly climatology of sea surface
C       temperature, e.g. the Reynolds climatology.
C
C   >>> ALLOW_CLIMSSS_RELAXATION <<<
C       Allow the relaxation to a monthly climatology of sea surface
C       salinity, e.g. the Levitus climatology.
C
C   >>> USE_EXF_INTERPOLATION <<<
C       Allows to provide input field on arbitrary Lat-Lon input grid
C       (as specified in EXF_NML_04) and to interpolate to model grid.
C     Note: default is to interpolate unless {FLD}_interpMethod is set to 0
C
C   ====================================================================
C
C    The following CPP options:
C       ALLOW_ATM_WIND / useAtmWind (useWind)
C       ALLOW_ATM_TEMP               (TEMP)
C       ALLOW_DOWNWARD_RADIATION     (DOWN)
C       ALLOW_BULKFORMULAE           (BULK)
C       EXF_READ_EVAP                (EVAP)
C       ALLOW_READ_TURBFLUXES        (TURB)
C
C    permit all ocean-model forcing configurations listed in the 2 tables below.
C    The first configuration (A1,B1) is the flux-forced, ocean model.
C    Configurations A2,B3 and A2,B4 use pkg/exf open-water bulk formulae
C    to compute, from atmospheric variables, the missing surface fluxes.
C    The forcing fields in the rightmost column are defined in EXF_FIELDS.h
C    (ocean-model surface forcing field are defined in model/inc/FFIELDS.h)
C
C    (A) Surface momentum flux: [model: fu,fv ; exf: ustress,vstress]
C
C    # |useWind|        actions
C   ---|-------|-------------------------------------------------------------
C   (1)| False | Read-in ustress,vstress (if needed in B, compute wind-speed)
C      |       |
C   (2)| True  | Read-in uwind,vwind ; compute wind stress ustress,vstress.
C   ---|-------|-------------------------------------------------------------
C
C    (B) Surface buoyancy flux:
C        [ net heat flux: Qnet (exf: hflux), net short-wave: Qsw (exf: swflux)
C          fresh-water flux: EmPmR (exf: sflux) and saltFlux (exf: saltflx) ]
C
C    # |TEMP |DOWN |BULK |EVAP |TURB |            actions
C   ---|-----|-----|-----|-----|-----|-------------------------------------
C   (1)|  -  |  -  |  -  |  -  |  -  | Read-in hflux, swflux and sflux.
C      |     |     |     |     |     |
C   (2)|  -  | def |  -  |  -  |  -  | Read-in hflux, swdown and sflux.
C      |     |     |     |     |     | Compute swflux.
C      |     |     |     |     |     |
C   (3)| def | def | def |  -  |  -  | Read-in atemp, aqh, swdown, lwdown,
C      |     |     |     |     |     |  precip, and runoff.
C      |     |     |     |     |     | Compute hflux, swflux and sflux.
C      |     |     |     |     |     |
C   (4)| def |  -  | def |  -  |  -  | Read-in atemp, aqh, swflux, lwflux,
C      |     |     |     |     |     |  precip, and runoff.
C      |     |     |     |     |     | Compute hflux and sflux.
C      |     |     |     |     |     |
C   (5)| def | def |  -  | def | def | Read-in hs, hl, swdown, lwdown,
C      |     |     |     |     |     |  evap, precip and runoff.
C      |     |     |     |     |     | Compute hflux, swflux and sflux.
C      |     |     |     |     |     |
C   (6)| def |  -  |  -  | def | def | Read-in hs, hl, swflux, lwflux,
C      |     |     |     |     |     |  evap, precip and runoff.
C      |     |     |     |     |     | Compute  hflux and sflux.
C
C   =======================================================================

C-  Bulk formulae related flags.
C Note: To use ALLOW_BULKFORMULAE or EXF_READ_EVAP, needs #define 

C-  Other forcing fields

C Note: To use EXF_CALC_ATMRHO, both 
C       and  need to be defined

C-  Zenith Angle/Albedo related flags.

C-  Use ocean_emissivity*lwdown in lwFlux. This flag should be defined
C   unless to reproduce old results (obtained with inconsistent old code)

C-  Relaxation to monthly climatologies.

C-  Allows to read-in (2-d) tidal geopotential forcing

C-  Allows to read-in seaice fraction from files (areaMaskFile)

C-  Use spatial interpolation to interpolate
C   forcing files from input grid to model grid.
C   for interpolated vector fields, rotate towards model-grid axis
C   using old rotation formulae (instead of grid-angles)
C   for interpolation around N & S pole, use the old formulation
C   (no pole symmetry, single vector-comp interp, reset to 0 zonal-comp @ N.pole)



CBOP
C     !ROUTINE: EXF_DIAGNOSTICS_INIT
C     !INTERFACE:
      SUBROUTINE EXF_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE EXF_DIAGNOSTICS_INIT
C     | o Routine to initialize list of all available diagnostics
C     |   for EXF package
C     *==========================================================*
C     \ev
C     !USES:
      IMPLICIT NONE

C     === Global variables ===
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
C
C     ==================================================================
C     HEADER EXF_PARAM.h
C     ==================================================================
C
C     o Header file for the surface flux data. Used by the external
C       forcing package.
C
C     started: Christian Eckert eckert@mit.edu  30-Jun-1999
C
C     changed: Christian Eckert eckert@mit.edu  14-Jan-2000
C              - Restructured the original version in order to have a
C                better interface to the MITgcmUV.
C
C              Christian Eckert eckert@mit.edu  12-Feb-2000
C              - Changed some variables names (package prefix: exf_)
C
C              Patrick Heimbach, heimbach@mit.edu  04-May-2000
C              - included exf_iprec to enable easy
C                switch between 32bit/64 bit data format
C
C              Patrick Heimbach, heimbach@mit.edu  01-May-2001
C              - added obcs parameters
C
C     mods for pkg/seaice: menemenlis@jpl.nasa.gov 20-Dec-2002
C
C     ==================================================================
C     HEADER EXF_PARAM.h
C     ==================================================================

C     Repeat period for forcing fields (s)
C     For example, for yearly repeat period: repeatPeriod=31556925.
C     Note: this option is not yet coded for sub-daily
C           forcing and for leap years but this limitation can be
C           circumvented by using a 4-year (1461-day) repeatPeriod
      Real*8     repeatPeriod

C     useExfCheckRange   :: check range of input/output field values
C     useExfYearlyFields :: when set, automatically add extension
C                           _YEAR to input file names; the yearly files need
C                           to contain all the records that pertain to
C                           a particular year, including day 1, hour zero
C     twoDigitYear       :: when set, use 2-digit year extension YR
C                           instead of _YEAR for useExfYearlyFields
C    useOBCSYearlyFields :: when reading Open-Boundary values, assume yearly
C                           climatology (def=false)
C     readStressOnAgrid  :: read wind-streess located on model-grid,
C                            A-grid position
C     rotateStressOnAgrid:: rotate from zonal/meridional components to
C                           U/V components
C     readStressOnCgrid  :: read wind-streess located on model-grid, C-grid
C                           position
C     stressIsOnCgrid    :: ustress & vstress are positioned on Arakawa C-grid
C     useAtmWind         :: use wind vector (uwind/vwind) to compute
C                           the wind stress (ustress/vstress)
C     useRelativeWind    :: Subtract U/VVEL or U/VICE from U/VWIND before
C                           computing U/VSTRESS
C     noNegativeEvap     :: prevent negative evap (= sea-surface condensation)
C     useStabilityFct_overIce :: over sea-ice, compute turbulent transfert
C                           coeff. function of stability (like over
C                           open ocean) rather than using fixed Coeff.
C     diags_opOceWeighted:: weight surface flux diagnostics with open-ocean
C                           fraction
C     useExfZenAlbedo    :: ocean albedo (direct part) may vary
C                           with zenith angle (see select_ZenAlbedo)
C     select_ZenAlbedo   :: switch to different methods to compute albedo
C                           (direct part)
C                        :: 0 just use exf_albedo
C                        :: 1 use daily mean albedo from exf_zenithangle_table.F
C                        :: 2 use daily mean albedo computed as in pkg/aim_v23
C                        :: 3 use daily variable albedo
C     useExfZenIncoming  :: compute incoming solar radiation along with
C                           zenith angle
C     exf_debugLev       :: select message printing to STDOUT (e.g., when
C                           read rec)
C     exf_monFreq        :: Monitor Frequency (s) for EXF
C     exf_adjMonFreq     :: Monitor Frequency (s) for AD exf variables
C     exf_adjMonSelect   :: select group of exf AD-variables to monitor
C                           =0 : none
C                           =1 : ocean forcing fu, fv, qnet, empmr (default)
C                           =2 : + atmospheric forcing fields (u/vwind,
C                                  atemp, lwdown, precip, etc.)
C                           =3 : + derived forcing fields (u/vstress,
C                                  h/sflux, wspeed)

      LOGICAL useExfCheckRange
      LOGICAL useExfYearlyFields, twoDigitYear
      LOGICAL useOBCSYearlyFields
      LOGICAL readStressOnAgrid
      LOGICAL rotateStressOnAgrid
      LOGICAL readStressOnCgrid
      LOGICAL stressIsOnCgrid
      LOGICAL useAtmWind
      LOGICAL useRelativeWind
      LOGICAL noNegativeEvap
      LOGICAL useStabilityFct_overIce
      LOGICAL diags_opOceWeighted

      LOGICAL useExfZenAlbedo
      INTEGER select_ZenAlbedo
      LOGICAL useExfZenIncoming

      INTEGER exf_debugLev
      Real*8     exf_monFreq
      Real*8     exf_adjMonFreq
      INTEGER exf_adjMonSelect

C     Drag coefficient scaling factor
      Real*8     exf_scal_BulkCdn

C     Maximum absolute windstress, used to reset unreastically high
C     data values
      Real*8     windstressmax

C     freezing temperature is the minimum temperature allowed, used
C     to reset climatological temperatures fields where they have
C     values below climtempfreeze
      Real*8 climtempfreeze

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C     Description of contents of surface boundary condition files
C     Note: fieldperiod=0 means input file is one time-constant field
C           fieldperiod=-12 means input file contains 12 monthly means
C-    for each field:
C     {fld}file       :: file-name for this field
C     {fld}startdate1 :: field starting date (YYYYMMDD)
C     {fld}startdate1 :: field starting date (YYYYMMDD)
C     {fld}startdate2 :: field starting date (HHMMSS)
C     {fld}StartTime  :: corresponding starting time (in sec) for this field
C     {fld}period     :: time period (in sec) between 2 reccords
C     {fld}RepCycle   :: time duration of a repeating cycle
C     {fld}const      :: uniform default field value

      INTEGER hfluxstartdate1
      INTEGER hfluxstartdate2
      Real*8     hfluxStartTime
      Real*8     hfluxperiod
      Real*8     hfluxRepCycle
      Real*8     hfluxconst
      Real*8     hflux_exfremo_intercept
      Real*8     hflux_exfremo_slope
      CHARACTER*1 hfluxmask

      INTEGER atempstartdate1
      INTEGER atempstartdate2
      Real*8     atempStartTime
      Real*8     atempperiod
      Real*8     atempRepCycle
      Real*8     atempconst
      Real*8     atemp_exfremo_intercept
      Real*8     atemp_exfremo_slope
      CHARACTER*1 atempmask

      INTEGER aqhstartdate1
      INTEGER aqhstartdate2
      Real*8     aqhStartTime
      Real*8     aqhperiod
      Real*8     aqhRepCycle
      Real*8     aqhconst
      Real*8     aqh_exfremo_intercept
      Real*8     aqh_exfremo_slope
      CHARACTER*1 aqhmask

      INTEGER hs_startdate1
      INTEGER hs_startdate2
      Real*8     hs_StartTime
      Real*8     hs_period
      Real*8     hs_RepCycle
      Real*8     hs_const
      Real*8     hs_exfremo_intercept
      Real*8     hs_exfremo_slope
      CHARACTER*1 hs_mask

      INTEGER hl_startdate1
      INTEGER hl_startdate2
      Real*8     hl_StartTime
      Real*8     hl_period
      Real*8     hl_RepCycle
      Real*8     hl_const
      Real*8     hl_exfremo_intercept
      Real*8     hl_exfremo_slope
      CHARACTER*1 hl_mask

      INTEGER sfluxstartdate1
      INTEGER sfluxstartdate2
      Real*8     sfluxStartTime
      Real*8     sfluxperiod
      Real*8     sfluxRepCycle
      Real*8     sfluxconst
      Real*8     sflux_exfremo_intercept
      Real*8     sflux_exfremo_slope
      CHARACTER*1 sfluxmask

      INTEGER evapstartdate1
      INTEGER evapstartdate2
      Real*8     evapStartTime
      Real*8     evapperiod
      Real*8     evapRepCycle
      Real*8     evapconst
      Real*8     evap_exfremo_intercept
      Real*8     evap_exfremo_slope
      CHARACTER*1 evapmask

      INTEGER precipstartdate1
      INTEGER precipstartdate2
      Real*8     precipStartTime
      Real*8     precipperiod
      Real*8     precipRepCycle
      Real*8     precipconst
      Real*8     precip_exfremo_intercept
      Real*8     precip_exfremo_slope
      CHARACTER*1 precipmask

      INTEGER snowprecipstartdate1
      INTEGER snowprecipstartdate2
      Real*8     snowprecipStartTime
      Real*8     snowprecipperiod
      Real*8     snowprecipRepCycle
      Real*8     snowprecipconst
      Real*8     snowprecip_exfremo_intercept
      Real*8     snowprecip_exfremo_slope
      CHARACTER*1 snowprecipmask

      INTEGER runoffstartdate1
      INTEGER runoffstartdate2
      Real*8     runoffStartTime
      Real*8     runoffperiod
      Real*8     runoffRepCycle
      Real*8     runoffconst
      Real*8     runoff_exfremo_intercept
      Real*8     runoff_exfremo_slope
      CHARACTER*1 runoffmask

      Real*8     runoftempconst
      Real*8     runoftemp_exfremo_intercept
      Real*8     runoftemp_exfremo_slope

      INTEGER saltflxstartdate1
      INTEGER saltflxstartdate2
      Real*8     saltflxStartTime
      Real*8     saltflxperiod
      Real*8     saltflxRepCycle
      Real*8     saltflxconst
      Real*8     saltflx_exfremo_intercept
      Real*8     saltflx_exfremo_slope
      CHARACTER*1 saltflxmask

      INTEGER ustressstartdate1
      INTEGER ustressstartdate2
      Real*8     ustressStartTime
      Real*8     ustressperiod
      Real*8     ustressRepCycle
      Real*8     ustressconst
      Real*8     ustress_exfremo_intercept
      Real*8     ustress_exfremo_slope
      CHARACTER*1 ustressmask

      INTEGER vstressstartdate1
      INTEGER vstressstartdate2
      Real*8     vstressStartTime
      Real*8     vstressperiod
      Real*8     vstressRepCycle
      Real*8     vstressconst
      Real*8     vstress_exfremo_intercept
      Real*8     vstress_exfremo_slope
      CHARACTER*1 vstressmask

      INTEGER uwindstartdate1
      INTEGER uwindstartdate2
      Real*8     uwindStartTime
      Real*8     uwindperiod
      Real*8     uwindRepCycle
      Real*8     uwindconst
      Real*8     uwind_exfremo_intercept
      Real*8     uwind_exfremo_slope
      CHARACTER*1 uwindmask

      INTEGER vwindstartdate1
      INTEGER vwindstartdate2
      Real*8     vwindStartTime
      Real*8     vwindperiod
      Real*8     vwindRepCycle
      Real*8     vwindconst
      Real*8     vwind_exfremo_intercept
      Real*8     vwind_exfremo_slope
      CHARACTER*1 vwindmask

      INTEGER wspeedstartdate1
      INTEGER wspeedstartdate2
      Real*8     wspeedStartTime
      Real*8     wspeedperiod
      Real*8     wspeedRepCycle
      Real*8     wspeedconst
      Real*8     wspeed_exfremo_intercept
      Real*8     wspeed_exfremo_slope
      CHARACTER*1 wspeedmask

      INTEGER swfluxstartdate1
      INTEGER swfluxstartdate2
      Real*8     swfluxStartTime
      Real*8     swfluxperiod
      Real*8     swfluxRepCycle
      Real*8     swfluxconst
      Real*8     swflux_exfremo_intercept
      Real*8     swflux_exfremo_slope
      CHARACTER*1 swfluxmask

      INTEGER lwfluxstartdate1
      INTEGER lwfluxstartdate2
      Real*8     lwfluxStartTime
      Real*8     lwfluxperiod
      Real*8     lwfluxRepCycle
      Real*8     lwfluxconst
      Real*8     lwflux_exfremo_intercept
      Real*8     lwflux_exfremo_slope
      CHARACTER*1 lwfluxmask

      INTEGER swdownstartdate1
      INTEGER swdownstartdate2
      Real*8     swdownStartTime
      Real*8     swdownperiod
      Real*8     swdownRepCycle
      Real*8     swdownconst
      Real*8     swdown_exfremo_intercept
      Real*8     swdown_exfremo_slope
      CHARACTER*1 swdownmask

      INTEGER lwdownstartdate1
      INTEGER lwdownstartdate2
      Real*8     lwdownStartTime
      Real*8     lwdownperiod
      Real*8     lwdownRepCycle
      Real*8     lwdownconst
      Real*8     lwdown_exfremo_intercept
      Real*8     lwdown_exfremo_slope
      CHARACTER*1 lwdownmask

      INTEGER apressurestartdate1
      INTEGER apressurestartdate2
      Real*8     apressureStartTime
      Real*8     apressureperiod
      Real*8     apressureRepCycle
      Real*8     apressureconst
      Real*8     apressure_exfremo_intercept
      Real*8     apressure_exfremo_slope
      CHARACTER*1 apressuremask

      INTEGER tidePotStartdate1
      INTEGER tidePotStartdate2
      Real*8     tidePotStartTime
      Real*8     tidePotPeriod
      Real*8     tidePotRepCycle
      Real*8     tidePotConst
      Real*8     tidePot_exfremo_intercept
      Real*8     tidePot_exfremo_slope
      CHARACTER*1 tidePotMask

      INTEGER areamaskstartdate1
      INTEGER areamaskstartdate2
      Real*8     areamaskStartTime
      Real*8     areamaskperiod
      Real*8     areamaskRepCycle
      Real*8     areamaskTauRelax
      Real*8     areamaskconst
      Real*8     areamask_exfremo_intercept
      Real*8     areamask_exfremo_slope
      CHARACTER*1 areamaskmask

C     Calendar data.
      INTEGER climsststartdate1
      INTEGER climsststartdate2
      Real*8     climsstStartTime
      Real*8     climsstperiod
      Real*8     climsstRepCycle
      Real*8     climsstTauRelax
      Real*8     climsstconst
      Real*8     climsst_exfremo_intercept
      Real*8     climsst_exfremo_slope
      CHARACTER*1 climsstmask

      INTEGER climsssstartdate1
      INTEGER climsssstartdate2
      Real*8     climsssStartTime
      Real*8     climsssperiod
      Real*8     climsssRepCycle
      Real*8     climsssTauRelax
      Real*8     climsssconst
      Real*8     climsss_exfremo_intercept
      Real*8     climsss_exfremo_slope
      CHARACTER*1 climsssmask

      INTEGER climustrstartdate1
      INTEGER climustrstartdate2
      Real*8     climustrStartTime
      Real*8     climustrperiod
      Real*8     climustrRepCycle
      Real*8     climustrTauRelax
      Real*8     climustrconst
      Real*8     climustr_exfremo_intercept
      Real*8     climustr_exfremo_slope
      CHARACTER*1 climustrmask

      INTEGER climvstrstartdate1
      INTEGER climvstrstartdate2
      Real*8     climvstrStartTime
      Real*8     climvstrperiod
      Real*8     climvstrRepCycle
      Real*8     climvstrTauRelax
      Real*8     climvstrconst
      Real*8     climvstr_exfremo_intercept
      Real*8     climvstr_exfremo_slope
      CHARACTER*1 climvstrmask

C-    The following variables are used in conjunction with pkg/obcs
C     to describe S/T/U/V open boundary condition files
      INTEGER obcsNstartdate1
      INTEGER obcsNstartdate2
      INTEGER obcsSstartdate1
      INTEGER obcsSstartdate2
      INTEGER obcsEstartdate1
      INTEGER obcsEstartdate2
      INTEGER obcsWstartdate1
      INTEGER obcsWstartdate2
      Real*8     obcsNstartTime
      Real*8     obcsNperiod
      Real*8     obcsNrepCycle
      Real*8     obcsSstartTime
      Real*8     obcsSperiod
      Real*8     obcsSrepCycle
      Real*8     obcsEstartTime
      Real*8     obcsEperiod
      Real*8     obcsErepCycle
      Real*8     obcsWstartTime
      Real*8     obcsWperiod
      Real*8     obcsWrepCycle

C-    The following variables are used in conjunction with pkg/obcs
C     and pkg/seaice to describe area, heff, hsnow, hsalt, uice,
C     and vice open boundary condition files
      INTEGER siobNstartdate1
      INTEGER siobNstartdate2
      INTEGER siobSstartdate1
      INTEGER siobSstartdate2
      INTEGER siobEstartdate1
      INTEGER siobEstartdate2
      INTEGER siobWstartdate1
      INTEGER siobWstartdate2
      Real*8     siobNstartTime
      Real*8     siobNperiod
      Real*8     siobNrepCycle
      Real*8     siobSstartTime
      Real*8     siobSperiod
      Real*8     siobSrepCycle
      Real*8     siobEstartTime
      Real*8     siobEperiod
      Real*8     siobErepCycle
      Real*8     siobWstartTime
      Real*8     siobWperiod
      Real*8     siobWrepCycle

C-    File names.
      CHARACTER*(128) hfluxfile
      CHARACTER*(128) atempfile
      CHARACTER*(128) aqhfile
      CHARACTER*(128) hs_file
      CHARACTER*(128) hl_file
      CHARACTER*(128) evapfile
      CHARACTER*(128) precipfile
      CHARACTER*(128) snowprecipfile
      CHARACTER*(128) sfluxfile
      CHARACTER*(128) runofffile
      CHARACTER*(128) runoftempfile
      CHARACTER*(128) saltflxfile
      CHARACTER*(128) ustressfile
      CHARACTER*(128) vstressfile
      CHARACTER*(128) uwindfile
      CHARACTER*(128) vwindfile
      CHARACTER*(128) wspeedfile
      CHARACTER*(128) swfluxfile
      CHARACTER*(128) lwfluxfile
      CHARACTER*(128) swdownfile
      CHARACTER*(128) lwdownfile
      CHARACTER*(128) apressurefile
      CHARACTER*(128) tidePotFile
      CHARACTER*(128) areamaskfile
      CHARACTER*(128) climsstfile
      CHARACTER*(128) climsssfile
      CHARACTER*(128) climustrfile
      CHARACTER*(128) climvstrfile

      COMMON /EXF_PARAM_L/
     &       useExfCheckRange,
     &       useExfYearlyFields, twoDigitYear,
     &       useOBCSYearlyFields,
     &       useExfZenAlbedo, useExfZenIncoming,
     &       readStressOnAgrid, readStressOnCgrid,
     &       stressIsOnCgrid, rotateStressOnAgrid,
     &       useAtmWind, useRelativeWind, noNegativeEvap,
     &       useStabilityFct_overIce, diags_opOceWeighted

      COMMON /EXF_PARAM_I/
     &       select_ZenAlbedo,  exf_debugLev,    exf_adjMonSelect,
     &       hfluxstartdate1,   hfluxstartdate2,
     &       atempstartdate1,   atempstartdate2,
     &       aqhstartdate1,     aqhstartdate2,
     &       hs_startdate1,     hs_startdate2,
     &       hl_startdate1,     hl_startdate2,
     &       sfluxstartdate1,   sfluxstartdate2,
     &       evapstartdate1,    evapstartdate2,
     &       runoffstartdate1,  runoffstartdate2,
     &       saltflxstartdate1, saltflxstartdate2,
     &       precipstartdate1,  precipstartdate2,
     &       snowprecipstartdate1, snowprecipstartdate2,
     &       ustressstartdate1, ustressstartdate2,
     &       vstressstartdate1, vstressstartdate2,
     &       uwindstartdate1,   uwindstartdate2,
     &       vwindstartdate1,   vwindstartdate2,
     &       wspeedstartdate1,  wspeedstartdate2,
     &       swfluxstartdate1,  swfluxstartdate2,
     &       lwfluxstartdate1,  lwfluxstartdate2,
     &       swdownstartdate1,  swdownstartdate2,
     &       lwdownstartdate1,  lwdownstartdate2,
     &       apressurestartdate1, apressurestartdate2,
     &       tidePotStartdate1, tidePotStartdate2,
     &       areamaskstartdate1,  areamaskstartdate2,
     &       obcsNstartdate1,   obcsNstartdate2,
     &       obcsSstartdate1,   obcsSstartdate2,
     &       obcsEstartdate1,   obcsEstartdate2,
     &       obcsWstartdate1,   obcsWstartdate2,
     &       siobNstartdate1,   siobNstartdate2,
     &       siobSstartdate1,   siobSstartdate2,
     &       siobEstartdate1,   siobEstartdate2,
     &       siobWstartdate1,   siobWstartdate2

      COMMON /EXF_PARAM_R/
     &       repeatPeriod,      exf_monFreq,     exf_adjMonFreq,
     &       exf_scal_BulkCdn,  windstressmax,
     &       hfluxconst,        hfluxRepCycle,
     &       hfluxperiod,       hfluxStartTime,
     &       atempconst,        atempRepCycle,
     &       atempperiod,       atempStartTime,
     &       aqhconst,          aqhRepCycle,
     &       aqhperiod,         aqhStartTime,
     &       hs_const,          hs_RepCycle,
     &       hs_period,         hs_StartTime,
     &       hl_const,          hl_RepCycle,
     &       hl_period,         hl_StartTime,
     &       sfluxconst,        sfluxRepCycle,
     &       sfluxperiod,       sfluxStartTime,
     &       evapconst,         evapRepCycle,
     &       evapperiod,        evapStartTime,
     &       precipconst,       precipRepCycle,
     &       precipperiod,      precipStartTime,
     &       snowprecipconst,   snowprecipRepCycle,
     &       snowprecipperiod,  snowprecipStartTime,
     &       runoffconst,       runoffRepCycle,
     &       runoffperiod,      runoffStartTime,
     &       runoftempconst,
     &       saltflxconst,      saltflxRepCycle,
     &       saltflxperiod,     saltflxStartTime,
     &       ustressconst,      ustressRepCycle,
     &       ustressperiod,     ustressStartTime,
     &       vstressconst,      vstressRepCycle,
     &       vstressperiod,     vstressStartTime,
     &       uwindconst,        uwindRepCycle,
     &       uwindperiod,       uwindStartTime,
     &       vwindconst,        vwindRepCycle,
     &       vwindperiod,       vwindStartTime,
     &       wspeedconst,       wspeedRepCycle,
     &       wspeedperiod,      wspeedStartTime,
     &       swfluxconst,       swfluxRepCycle,
     &       swfluxperiod,      swfluxStartTime,
     &       lwfluxconst,       lwfluxRepCycle,
     &       lwfluxperiod,      lwfluxStartTime,
     &       swdownconst,       swdownRepCycle,
     &       swdownperiod,      swdownStartTime,
     &       lwdownconst,       lwdownRepCycle,
     &       lwdownperiod,      lwdownStartTime,
     &       apressureconst,    apressureRepCycle,
     &       apressureperiod,   apressureStartTime,
     &       tidePotConst,      tidePotRepCycle,
     &       tidePotPeriod,     tidePotStartTime,
     &       areamaskconst,     areamaskRepCycle,
     &       areamaskperiod,    areamaskStartTime,
     &       obcsNrepCycle,     obcsNperiod,     obcsNstartTime,
     &       obcsSrepCycle,     obcsSperiod,     obcsSstartTime,
     &       obcsErepCycle,     obcsEperiod,     obcsEstartTime,
     &       obcsWrepCycle,     obcsWperiod,     obcsWstartTime,
     &       siobNrepCycle,     siobNperiod,     siobNstartTime,
     &       siobSrepCycle,     siobSperiod,     siobSstartTime,
     &       siobErepCycle,     siobEperiod,     siobEstartTime,
     &       siobWrepCycle,     siobWperiod,     siobWstartTime

      COMMON /EXF_PARAM_TREND_REMOVAL/
     &       hflux_exfremo_intercept,
     &       atemp_exfremo_intercept,
     &       aqh_exfremo_intercept,
     &       hs_exfremo_intercept,
     &       hl_exfremo_intercept,
     &       sflux_exfremo_intercept,
     &       evap_exfremo_intercept,
     &       precip_exfremo_intercept,
     &       snowprecip_exfremo_intercept,
     &       runoff_exfremo_intercept,
     &       runoftemp_exfremo_intercept,
     &       saltflx_exfremo_intercept,
     &       ustress_exfremo_intercept,
     &       vstress_exfremo_intercept,
     &       uwind_exfremo_intercept,
     &       vwind_exfremo_intercept,
     &       wspeed_exfremo_intercept,
     &       swflux_exfremo_intercept,
     &       lwflux_exfremo_intercept,
     &       swdown_exfremo_intercept,
     &       lwdown_exfremo_intercept,
     &       apressure_exfremo_intercept,
     &       tidePot_exfremo_intercept,
     &       areamask_exfremo_intercept,
     &       hflux_exfremo_slope,
     &       atemp_exfremo_slope,
     &       aqh_exfremo_slope,
     &       hs_exfremo_slope,
     &       hl_exfremo_slope,
     &       sflux_exfremo_slope,
     &       evap_exfremo_slope,
     &       precip_exfremo_slope,
     &       snowprecip_exfremo_slope,
     &       runoff_exfremo_slope,
     &       runoftemp_exfremo_slope,
     &       saltflx_exfremo_slope,
     &       ustress_exfremo_slope,
     &       vstress_exfremo_slope,
     &       uwind_exfremo_slope,
     &       vwind_exfremo_slope,
     &       wspeed_exfremo_slope,
     &       swflux_exfremo_slope,
     &       lwflux_exfremo_slope,
     &       swdown_exfremo_slope,
     &       lwdown_exfremo_slope,
     &       apressure_exfremo_slope,
     &       tidePot_exfremo_slope,
     &       areamask_exfremo_slope

      COMMON /EXF_PARAM_C/
     &       hfluxfile,     hfluxmask,
     &       atempfile,     atempmask,
     &       aqhfile,       aqhmask,
     &       hs_file,       hs_mask,
     &       hl_file,       hl_mask,
     &       sfluxfile,     sfluxmask,
     &       evapfile,      evapmask,
     &       precipfile,    precipmask,
     &       snowprecipfile,snowprecipmask,
     &       runofffile,    runoffmask,
     &       runoftempfile,
     &       saltflxfile,   saltflxmask,
     &       ustressfile,   ustressmask,
     &       vstressfile,   vstressmask,
     &       uwindfile,     uwindmask,
     &       vwindfile,     vwindmask,
     &       wspeedfile,    wspeedmask,
     &       swfluxfile,    swfluxmask,
     &       lwfluxfile,    lwfluxmask,
     &       swdownfile,    swdownmask,
     &       lwdownfile,    lwdownmask,
     &       apressurefile, apressuremask,
     &       tidePotFile,   tidePotMask,
     &       areamaskfile,  areamaskmask

      COMMON /EXF_CLIM_I/
     &       climsststartdate1,  climsststartdate2,
     &       climsssstartdate1,  climsssstartdate2,
     &       climustrstartdate1,  climustrstartdate2,
     &       climvstrstartdate1,  climvstrstartdate2

      COMMON /EXF_CLIM_C/
     &       climsstfile,  climsstmask,
     &       climsssfile,  climsssmask,
     &       climustrfile, climustrmask,
     &       climvstrfile, climvstrmask

      COMMON /EXF_CLIM_R/
     &       climtempfreeze,
     &       climsstconst,       climsstRepCycle,
     &       climsstperiod,      climsstStartTime,
     &       climsssconst,       climsssRepCycle,
     &       climsssperiod,      climsssStartTime,
     &       climustrconst,      climustrRepCycle,
     &       climustrperiod,     climustrStartTime,
     &       climvstrconst,      climvstrRepCycle,
     &       climvstrperiod,     climvstrStartTime,
     &       climsstTauRelax,    climsssTauRelax,
     &       climustrTauRelax,   climvstrTauRelax,
     &       areamaskTauRelax,
     &       climsst_exfremo_intercept, climsst_exfremo_slope,
     &       climsss_exfremo_intercept, climsss_exfremo_slope,
     &       climustr_exfremo_intercept, climustr_exfremo_slope,
     &       climvstr_exfremo_intercept, climvstr_exfremo_slope,
     &       exf_inscal_climsst, exf_inscal_climsss,
     &       exf_inscal_climustr, exf_inscal_climvstr

C     file precision and field type

      COMMON /EXF_PARAM_TYPE/
     &       exf_iprec,
     &       exf_iprec_obcs

      INTEGER exf_iprec
      INTEGER exf_iprec_obcs

C-    Scaling factors:
C     exf_inscal_{fld}   :: input scaling factors
C     exf_offset_atemp   :: input air temperature offset
C                        :: (for conversion from C to K, if needed)
C     exf_outscale_{fld} :: output scaling factors

      Real*8     exf_inscal_hflux
      Real*8     exf_inscal_sflux
      Real*8     exf_inscal_ustress
      Real*8     exf_inscal_vstress
      Real*8     exf_inscal_uwind
      Real*8     exf_inscal_vwind
      Real*8     exf_inscal_wspeed
      Real*8     exf_inscal_swflux
      Real*8     exf_inscal_lwflux
      Real*8     exf_inscal_precip
      Real*8     exf_inscal_snowprecip
c     Real*8     exf_inscal_sst
c     Real*8     exf_inscal_sss
      Real*8     exf_inscal_atemp, exf_offset_atemp
      Real*8     exf_inscal_aqh
      Real*8     exf_inscal_hs
      Real*8     exf_inscal_hl
      Real*8     exf_inscal_evap
      Real*8     exf_inscal_apressure
      Real*8     exf_inscal_runoff
      Real*8     exf_inscal_runoftemp
      Real*8     exf_inscal_saltflx
      Real*8     exf_inscal_swdown
      Real*8     exf_inscal_lwdown
      Real*8     exf_inscal_tidePot
      Real*8     exf_inscal_areamask
      Real*8     exf_inscal_climsst
      Real*8     exf_inscal_climsss
      Real*8     exf_inscal_climustr
      Real*8     exf_inscal_climvstr

      Real*8     exf_outscal_hflux
      Real*8     exf_outscal_sflux
      Real*8     exf_outscal_ustress
      Real*8     exf_outscal_vstress
      Real*8     exf_outscal_swflux
      Real*8     exf_outscal_sst
      Real*8     exf_outscal_sss
      Real*8     exf_outscal_apressure
      Real*8     exf_outscal_tidePot
      Real*8     exf_outscal_areamask

      COMMON /EXF_PARAM_SCAL/
     &                      exf_inscal_hflux,
     &                      exf_inscal_sflux,
     &                      exf_inscal_ustress,
     &                      exf_inscal_vstress,
     &                      exf_inscal_uwind,
     &                      exf_inscal_vwind,
     &                      exf_inscal_wspeed,
     &                      exf_inscal_swflux,
     &                      exf_inscal_lwflux,
     &                      exf_inscal_precip,
     &                      exf_inscal_snowprecip,
c    &                      exf_inscal_sst,
c    &                      exf_inscal_sss,
     &                      exf_inscal_atemp, exf_offset_atemp,
     &                      exf_inscal_aqh,
     &                      exf_inscal_hs,
     &                      exf_inscal_hl,
     &                      exf_inscal_evap,
     &                      exf_inscal_apressure,
     &                      exf_inscal_runoff,
     &                      exf_inscal_runoftemp,
     &                      exf_inscal_saltflx,
     &                      exf_inscal_swdown,
     &                      exf_inscal_lwdown,
     &                      exf_inscal_tidePot,
     &                      exf_inscal_areamask,
     &                      exf_outscal_hflux,
     &                      exf_outscal_sflux,
     &                      exf_outscal_ustress,
     &                      exf_outscal_vstress,
     &                      exf_outscal_swflux,
     &                      exf_outscal_sst,
     &                      exf_outscal_sss,
     &                      exf_outscal_apressure,
     &                      exf_outscal_tidePot,
     &                      exf_outscal_areamask

C- note: pkg/exf Interpolation parameters (#ifdef USE_EXF_INTERPOLATION )
C   have been moved to specific header file: EXF_INTERP_PARAM.h

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     === Local variables ===
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C              cf head of S/R DIAGNOSTICS_INIT_EARLY or DIAGNOSTICS_MAIN_INIT
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER       diagNum
      INTEGER       diagMate
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      diagName  = 'EXFhs   '
      diagTitle = 'Sensible heat flux into ocean, >0 increases theta'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFhl   '
      diagTitle = 'Latent heat flux into ocean, >0 increases theta'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFlwnet'
      diagTitle = 'Net upward longwave radiation, >0 decreases theta'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFswnet'
      diagTitle = 'Net upward shortwave radiation, >0 decreases theta'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFlwdn '
      diagTitle = 'Downward longwave radiation, >0 increases theta'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFswdn '
      diagTitle = 'Downward shortwave radiation, >0 increases theta'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFqnet '
      diagTitle = 'Net upward heat flux (turb+rad), >0 decreases theta'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C--   surface fluxes:
      diagName  = 'EXFtaux '
      diagTitle = 'zonal surface wind stress, >0 increases uVel    '
      diagUnits = 'N/m^2           '
      IF ( stressIsOnCgrid ) THEN
       diagCode = 'UU      U1      '
      ELSE
       diagCode = 'UM      U1      '
      ENDIF
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFtauy '
      diagTitle = 'meridional surface wind stress, >0 increases vVel'
      diagUnits = 'N/m^2           '
      IF ( stressIsOnCgrid ) THEN
       diagCode = 'VV      U1      '
      ELSE
       diagCode = 'VM      U1      '
      ENDIF
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFuwind'
      diagTitle = 'zonal 10-m wind speed, >0 eastward'
      diagUnits = 'm/s             '
      diagCode  = 'UM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFvwind'
      diagTitle = 'meridional 10-m wind speed, >0 northward'
      diagUnits = 'm/s             '
      diagCode  = 'VM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFwspee'
      diagTitle = '10-m wind speed modulus ( >= 0 )'
      diagUnits = 'm/s             '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFatemp'
      diagTitle = 'surface (2-m) air temperature'
      diagUnits = 'degK            '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFaqh  '
      diagTitle = 'surface (2-m) specific humidity'
      diagUnits = 'kg/kg           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFevap '
      diagTitle = 'evaporation, > 0 increases salinity'
      diagUnits = 'm/s             '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFpreci'
      diagTitle = 'precipitation, > 0 decreases salinity'
      diagUnits = 'm/s             '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFsnow '
      diagTitle = 'snow precipitation, > 0 decreases salinity'
      diagUnits = 'm/s             '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFempmr'
      diagTitle = 'net upward freshwater flux, > 0 increases salinity'
      diagUnits = 'm/s             '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFpress'
      diagTitle = 'atmospheric pressure field'
      diagUnits = 'N/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFroff '
      diagTitle = 'river runoff, > 0 decreases salinity'
      diagUnits = 'm/s             '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFroft '
      diagTitle = 'river runoff temperature'
      diagUnits = 'deg C           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EXFsalfx'
      diagTitle = 'upward salt flux, > 0 decreases salinity'
      diagUnits = 'g/m^2/s         '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )


      RETURN
      END
