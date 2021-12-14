












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



CBOP 0
C     !ROUTINE: EXF_WEIGHT_SFX_DIAGS
C     !INTERFACE:
      SUBROUTINE EXF_WEIGHT_SFX_DIAGS(
     I                                 siceFrac, diagGroup,
     I                                 myTime, myIter, myThid )

C     !DESCRIPTION:
C     Apply open-ocean fraction weight to surface fluxes
C     before filling-in corresponding diagnostics

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

C

CBOP
C     !ROUTINE: PARAMS.h
C     !INTERFACE:
C     #include PARAMS.h

C     !DESCRIPTION:
C     Header file defining model "parameters".  The values from the
C     model standard input file are stored into the variables held
C     here. Notes describing the parameters can also be found here.

CEOP

C--   Contants
C     Useful physical values
      Real*8 PI
      PARAMETER ( PI    = 3.14159265358979323844d0   )
      Real*8 deg2rad
      PARAMETER ( deg2rad = 2.d0*PI/360.d0           )

C--   COMMON /PARM_C/ Character valued parameters used by the model.
C     buoyancyRelation :: Flag used to indicate which relation to use to
C                         get buoyancy.
C     eosType         :: choose the equation of state:
C                        LINEAR, POLY3, UNESCO, JMD95Z, JMD95P, MDJWF, IDEALGAS
C     pickupSuff      :: force to start from pickup files (even if nIter0=0)
C                        and read pickup files with this suffix (max 10 Char.)
C     mdsioLocalDir   :: read-write tiled file from/to this directory name
C                        (+ 4 digits Processor-Rank) instead of current dir.
C     adTapeDir       :: read-write checkpointing tape files from/to this
C                        directory name instead of current dir. Conflicts
C                        mdsioLocalDir, so only one of the two can be set.
C     tRefFile      :: File containing reference Potential Temperat.  tRef (1.D)
C     sRefFile      :: File containing reference salinity/spec.humid. sRef (1.D)
C     rhoRefFile    :: File containing reference density profile rhoRef (1.D)
C     gravityFile   :: File containing gravity vertical profile (1.D)
C     delRFile      :: File containing vertical grid spacing delR  (1.D array)
C     delRcFile     :: File containing vertical grid spacing delRc (1.D array)
C     hybSigmFile   :: File containing hybrid-sigma vertical coord. coeff. (2x 1.D)
C     delXFile      :: File containing X-spacing grid definition (1.D array)
C     delYFile      :: File containing Y-spacing grid definition (1.D array)
C     horizGridFile :: File containing horizontal-grid definition
C                        (only when using curvilinear_grid)
C     bathyFile       :: File containing bathymetry. If not defined bathymetry
C                        is taken from inline function.
C     topoFile        :: File containing the topography of the surface (unit=m)
C                        (mainly used for the atmosphere = ground height).
C     addWwallFile    :: File containing 2-D additional Western  cell-edge wall
C     addSwallFile    :: File containing 2-D additional Southern cell-edge wall
C                        (e.g., to add "thin-wall" where it is =1)
C     hydrogThetaFile :: File containing initial hydrographic data (3-D)
C                        for potential temperature.
C     hydrogSaltFile  :: File containing initial hydrographic data (3-D)
C                        for salinity.
C     diffKrFile      :: File containing 3D specification of vertical diffusivity
C     viscAhDfile     :: File containing 3D specification of horizontal viscosity
C     viscAhZfile     :: File containing 3D specification of horizontal viscosity
C     viscA4Dfile     :: File containing 3D specification of horizontal viscosity
C     viscA4Zfile     :: File containing 3D specification of horizontal viscosity
C     zonalWindFile   :: File containing zonal wind data
C     meridWindFile   :: File containing meridional wind data
C     thetaClimFile   :: File containing surface theta climataology used
C                        in relaxation term -lambda(theta-theta*)
C     saltClimFile    :: File containing surface salt climataology used
C                        in relaxation term -lambda(salt-salt*)
C     surfQfile       :: File containing surface heat flux, excluding SW
C                        (old version, kept for backward compatibility)
C     surfQnetFile    :: File containing surface net heat flux
C     surfQswFile     :: File containing surface shortwave radiation
C     EmPmRfile       :: File containing surface fresh water flux
C           NOTE: for backward compatibility EmPmRfile is specified in
C                 m/s when using external_fields_load.F.  It is converted
C                 to kg/m2/s by multiplying by rhoConstFresh.
C     saltFluxFile    :: File containing surface salt flux
C     pLoadFile       :: File containing pressure loading
C     geoPotAnomFile  :: File containing constant geopotential anomaly due to
C                        density structure
C     addMassFile     :: File containing source/sink of fluid in the interior
C     eddyPsiXFile    :: File containing zonal Eddy streamfunction data
C     eddyPsiYFile    :: File containing meridional Eddy streamfunction data
C     the_run_name    :: string identifying the name of the model "run"
      COMMON /PARM_C/
     &                buoyancyRelation, eosType,
     &                pickupSuff, mdsioLocalDir, adTapeDir,
     &                tRefFile, sRefFile, rhoRefFile, gravityFile,
     &                delRFile, delRcFile, hybSigmFile,
     &                delXFile, delYFile, horizGridFile,
     &                bathyFile, topoFile, addWwallFile, addSwallFile,
     &                viscAhDfile, viscAhZfile,
     &                viscA4Dfile, viscA4Zfile,
     &                hydrogThetaFile, hydrogSaltFile, diffKrFile,
     &                zonalWindFile, meridWindFile, thetaClimFile,
     &                saltClimFile,
     &                EmPmRfile, saltFluxFile,
     &                surfQfile, surfQnetFile, surfQswFile,
     &                lambdaThetaFile, lambdaSaltFile,
     &                uVelInitFile, vVelInitFile, pSurfInitFile,
     &                pLoadFile, geoPotAnomFile, addMassFile,
     &                eddyPsiXFile, eddyPsiYFile, geothermalFile,
     &                the_run_name
      CHARACTER*(MAX_LEN_FNAM) buoyancyRelation
      CHARACTER*(6)  eosType
      CHARACTER*(10) pickupSuff
      CHARACTER*(MAX_LEN_FNAM) mdsioLocalDir
      CHARACTER*(MAX_LEN_FNAM) adTapeDir
      CHARACTER*(MAX_LEN_FNAM) tRefFile
      CHARACTER*(MAX_LEN_FNAM) sRefFile
      CHARACTER*(MAX_LEN_FNAM) rhoRefFile
      CHARACTER*(MAX_LEN_FNAM) gravityFile
      CHARACTER*(MAX_LEN_FNAM) delRFile
      CHARACTER*(MAX_LEN_FNAM) delRcFile
      CHARACTER*(MAX_LEN_FNAM) hybSigmFile
      CHARACTER*(MAX_LEN_FNAM) delXFile
      CHARACTER*(MAX_LEN_FNAM) delYFile
      CHARACTER*(MAX_LEN_FNAM) horizGridFile
      CHARACTER*(MAX_LEN_FNAM) bathyFile, topoFile
      CHARACTER*(MAX_LEN_FNAM) addWwallFile, addSwallFile
      CHARACTER*(MAX_LEN_FNAM) hydrogThetaFile, hydrogSaltFile
      CHARACTER*(MAX_LEN_FNAM) diffKrFile
      CHARACTER*(MAX_LEN_FNAM) viscAhDfile
      CHARACTER*(MAX_LEN_FNAM) viscAhZfile
      CHARACTER*(MAX_LEN_FNAM) viscA4Dfile
      CHARACTER*(MAX_LEN_FNAM) viscA4Zfile
      CHARACTER*(MAX_LEN_FNAM) zonalWindFile
      CHARACTER*(MAX_LEN_FNAM) meridWindFile
      CHARACTER*(MAX_LEN_FNAM) thetaClimFile
      CHARACTER*(MAX_LEN_FNAM) saltClimFile
      CHARACTER*(MAX_LEN_FNAM) surfQfile
      CHARACTER*(MAX_LEN_FNAM) surfQnetFile
      CHARACTER*(MAX_LEN_FNAM) surfQswFile
      CHARACTER*(MAX_LEN_FNAM) EmPmRfile
      CHARACTER*(MAX_LEN_FNAM) saltFluxFile
      CHARACTER*(MAX_LEN_FNAM) uVelInitFile
      CHARACTER*(MAX_LEN_FNAM) vVelInitFile
      CHARACTER*(MAX_LEN_FNAM) pSurfInitFile
      CHARACTER*(MAX_LEN_FNAM) pLoadFile
      CHARACTER*(MAX_LEN_FNAM) geoPotAnomFile
      CHARACTER*(MAX_LEN_FNAM) addMassFile
      CHARACTER*(MAX_LEN_FNAM) eddyPsiXFile
      CHARACTER*(MAX_LEN_FNAM) eddyPsiYFile
      CHARACTER*(MAX_LEN_FNAM) geothermalFile
      CHARACTER*(MAX_LEN_FNAM) lambdaThetaFile
      CHARACTER*(MAX_LEN_FNAM) lambdaSaltFile
      CHARACTER*(MAX_LEN_PREC/2) the_run_name

C--   COMMON /PARM_I/ Integer valued parameters used by the model.
C     cg2dMaxIters        :: Maximum number of iterations in the
C                            two-dimensional con. grad solver.
C     cg2dMinItersNSA     :: Minimum number of iterations in the
C                            not-self-adjoint version (cg2d_nsa.F) of the
C                            two-dimensional con. grad solver (default = 0).
C     cg2dPreCondFreq     :: Frequency for updating cg2d preconditioner
C                            (non-linear free-surf.)
C     cg2dUseMinResSol    :: =0 : use last-iteration/converged solution
C                            =1 : use solver minimum-residual solution
C     cg3dMaxIters        :: Maximum number of iterations in the
C                            three-dimensional con. grad solver.
C     printResidualFreq   :: Frequency for printing residual in CG iterations
C     nIter0              :: Start time-step number of for this run
C     nTimeSteps          :: Number of timesteps to execute
C     nTimeSteps_l2       :: Number of inner timesteps to execute per timestep
C     selectCoriMap       :: select setting of Coriolis parameter map:
C                           =0 f-Plane (Constant Coriolis, = f0)
C                           =1 Beta-Plane Coriolis (= f0 + beta.y)
C                           =2 Spherical Coriolis (= 2.omega.sin(phi))
C                           =3 Read Coriolis 2-d fields from files.
C     selectSigmaCoord    :: option related to sigma vertical coordinate
C     nonlinFreeSurf      :: option related to non-linear free surface
C                           =0 Linear free surface ; >0 Non-linear
C     select_rStar        :: option related to r* vertical coordinate
C                           =0 (default) use r coord. ; > 0 use r*
C     selectNHfreeSurf    :: option for Non-Hydrostatic (free-)Surface formulation:
C                           =0 (default) hydrostatic surf. ; > 0 add NH effects.
C     selectP_inEOS_Zc    :: select which pressure to use in EOS (for z-coords)
C                           =0: simply: -g*rhoConst*z
C                           =1: use pRef = integral{-g*rho(Tref,Sref,pRef)*dz}
C                           =2: use hydrostatic dynamical pressure
C                           =3: use full (Hyd+NH) dynamical pressure
C     selectAddFluid      :: option to add mass source/sink of fluid in the interior
C                            (3-D generalisation of oceanic real-fresh water flux)
C                           =0 off ; =1 add fluid ; =-1 virtual flux (no mass added)
C     selectImplicitDrag  :: select Implicit treatment of bottom/top drag
C                           = 0: fully explicit
C                           = 1: implicit on provisional velocity
C                                (i.e., before grad.Eta increment)
C                           = 2: fully implicit (combined with Impl Surf.Press)
C     momForcingOutAB     :: =1: take momentum forcing contribution
C                            out of (=0: in) Adams-Bashforth time stepping.
C     tracForcingOutAB    :: =1: take tracer (Temp,Salt,pTracers) forcing contribution
C                            out of (=0: in) Adams-Bashforth time stepping.
C     tempAdvScheme       :: Temp. Horiz.Advection scheme selector
C     tempVertAdvScheme   :: Temp. Vert. Advection scheme selector
C     saltAdvScheme       :: Salt. Horiz.advection scheme selector
C     saltVertAdvScheme   :: Salt. Vert. Advection scheme selector
C     selectKEscheme      :: Kinetic Energy scheme selector (Vector Inv.)
C     selectVortScheme    :: Scheme selector for Vorticity term (Vector Inv.)
C     selectCoriScheme    :: Scheme selector for Coriolis term
C     selectBotDragQuadr  :: quadratic bottom drag discretisation option:
C                           =0: average KE from grid center to U & V location
C                           =1: use local velocity norm @ U & V location
C                           =2: same with wet-point averaging of other component
C     pCellMix_select     :: select option to enhance mixing near surface & bottom
C                            unit digit: near bottom ; tens digit: near surface
C                            with digit =0 : disable ;
C                           = 1 : increases mixing linearly with recip_hFac
C                           = 2,3,4 : increases mixing by recip_hFac^(2,3,4)
C     readBinaryPrec      :: Precision used for reading binary files
C     writeStatePrec      :: Precision used for writing model state.
C     writeBinaryPrec     :: Precision used for writing binary files
C     rwSuffixType        :: controls the format of the mds file suffix.
C                          =0 (default): use iteration number (myIter, I10.10);
C                          =1: 100*myTime (100th sec); =2: myTime (seconds);
C                          =3: myTime/360 (10th of hr); =4: myTime/3600 (hours).
C     monitorSelect       :: select group of variables to monitor
C                            =1 : dynvars ; =2 : + vort ; =3 : + surface
C-    debugLevel          :: controls printing of algorithm intermediate results
C                            and statistics ; higher -> more writing
C-    plotLevel           :: controls printing of field maps ; higher -> more flds

      COMMON /PARM_I/
     &        cg2dMaxIters, cg2dMinItersNSA,
     &        cg2dPreCondFreq, cg2dUseMinResSol,
     &        cg3dMaxIters, printResidualFreq,
     &        nIter0, nTimeSteps, nTimeSteps_l2, nEndIter,
     &        selectCoriMap,
     &        selectSigmaCoord,
     &        nonlinFreeSurf, select_rStar,
     &        selectNHfreeSurf, selectP_inEOS_Zc,
     &        selectAddFluid, selectImplicitDrag,
     &        momForcingOutAB, tracForcingOutAB,
     &        tempAdvScheme, tempVertAdvScheme,
     &        saltAdvScheme, saltVertAdvScheme,
     &        selectKEscheme, selectVortScheme, selectCoriScheme,
     &        selectBotDragQuadr, pCellMix_select,
     &        readBinaryPrec, writeBinaryPrec, writeStatePrec,
     &        rwSuffixType, monitorSelect, debugLevel, plotLevel
      INTEGER cg2dMaxIters
      INTEGER cg2dMinItersNSA
      INTEGER cg2dPreCondFreq
      INTEGER cg2dUseMinResSol
      INTEGER cg3dMaxIters
      INTEGER printResidualFreq
      INTEGER nIter0
      INTEGER nTimeSteps
      INTEGER nTimeSteps_l2
      INTEGER nEndIter
      INTEGER selectCoriMap
      INTEGER selectSigmaCoord
      INTEGER nonlinFreeSurf
      INTEGER select_rStar
      INTEGER selectNHfreeSurf
      INTEGER selectP_inEOS_Zc
      INTEGER selectAddFluid
      INTEGER selectImplicitDrag
      INTEGER momForcingOutAB, tracForcingOutAB
      INTEGER tempAdvScheme, tempVertAdvScheme
      INTEGER saltAdvScheme, saltVertAdvScheme
      INTEGER selectKEscheme
      INTEGER selectVortScheme
      INTEGER selectCoriScheme
      INTEGER selectBotDragQuadr
      INTEGER pCellMix_select
      INTEGER readBinaryPrec
      INTEGER writeStatePrec
      INTEGER writeBinaryPrec
      INTEGER rwSuffixType
      INTEGER monitorSelect
      INTEGER debugLevel
      INTEGER plotLevel

C--   COMMON /PARM_L/ Logical valued parameters used by the model.
C- Coordinate + Grid params:
C     fluidIsAir       :: Set to indicate that the fluid major constituent
C                         is air
C     fluidIsWater     :: Set to indicate that the fluid major constituent
C                         is water
C     usingPCoords     :: Set to indicate that we are working in a pressure
C                         type coordinate (p or p*).
C     usingZCoords     :: Set to indicate that we are working in a height
C                         type coordinate (z or z*)
C     usingCartesianGrid :: If TRUE grid generation will be in a cartesian
C                           coordinate frame.
C     usingSphericalPolarGrid :: If TRUE grid generation will be in a
C                                spherical polar frame.
C     rotateGrid      :: rotate grid coordinates to geographical coordinates
C                        according to Euler angles phiEuler, thetaEuler, psiEuler
C     usingCylindricalGrid :: If TRUE grid generation will be Cylindrical
C     usingCurvilinearGrid :: If TRUE, use a curvilinear grid (to be provided)
C     hasWetCSCorners :: domain contains CS-type corners where dynamics is solved
C     deepAtmosphere :: deep model (drop the shallow-atmosphere approximation)
C     setInterFDr    :: set Interface depth (put cell-Center at the middle)
C     setCenterDr    :: set cell-Center depth (put Interface at the middle)
C     useMin4hFacEdges :: set hFacW,hFacS as minimum of adjacent hFacC factor
C     interViscAr_pCell :: account for partial-cell in interior vert. viscosity
C     interDiffKr_pCell :: account for partial-cell in interior vert. diffusion
C- Momentum params:
C     no_slip_sides  :: Impose "no-slip" at lateral boundaries.
C     no_slip_bottom :: Impose "no-slip" at bottom boundary.
C     bottomVisc_pCell :: account for partial-cell in bottom visc. (no-slip BC)
C     useSmag3D      :: Use isotropic 3-D Smagorinsky
C     useFullLeith   :: Set to true to use full Leith viscosity(may be unstable
C                       on irregular grids)
C     useStrainTensionVisc:: Set to true to use Strain-Tension viscous terms
C     useAreaViscLength :: Set to true to use old scaling for viscous lengths,
C                          e.g., L2=Raz.  May be preferable for cube sphere.
C     momViscosity  :: Flag which turns momentum friction terms on and off.
C     momAdvection  :: Flag which turns advection of momentum on and off.
C     momForcing    :: Flag which turns external forcing of momentum on and off.
C     momTidalForcing    :: Flag which turns tidal forcing on and off.
C     momPressureForcing :: Flag which turns pressure term in momentum equation
C                          on and off.
C     metricTerms   :: Flag which turns metric terms on or off.
C     useNHMTerms   :: If TRUE use non-hydrostatic metric terms.
C     useCoriolis   :: Flag which turns the coriolis terms on and off.
C     use3dCoriolis :: Turns the 3-D coriolis terms (in Omega.cos Phi) on - off
C     useCDscheme   :: use CD-scheme to calculate Coriolis terms.
C     vectorInvariantMomentum :: use Vector-Invariant form (mom_vecinv package)
C                                (default = F = use mom_fluxform package)
C     useJamartMomAdv :: Use wet-point method for V.I. non-linear term
C     upwindVorticity :: bias interpolation of vorticity in the Coriolis term
C     highOrderVorticity :: use 3rd/4th order interp. of vorticity (V.I., advection)
C     useAbsVorticity :: work with f+zeta in Coriolis terms
C     upwindShear     :: use 1rst order upwind interp. (V.I., vertical advection)
C     momStepping    :: Turns momentum equation time-stepping off
C     calc_wVelocity :: Turns vertical velocity calculation off
C- Temp. & Salt params:
C     tempStepping   :: Turns temperature equation time-stepping on/off
C     saltStepping   :: Turns salinity equation time-stepping on/off
C     addFrictionHeating :: account for frictional heating
C     temp_stayPositive :: use Smolarkiewicz Hack to ensure Temp stays positive
C     salt_stayPositive :: use Smolarkiewicz Hack to ensure Salt stays positive
C     tempAdvection  :: Flag which turns advection of temperature on and off.
C     tempVertDiff4  :: use vertical bi-harmonic diffusion for temperature
C     tempIsActiveTr :: Pot.Temp. is a dynamically active tracer
C     tempForcing    :: Flag which turns external forcing of temperature on/off
C     saltAdvection  :: Flag which turns advection of salinity on and off.
C     saltVertDiff4  :: use vertical bi-harmonic diffusion for salinity
C     saltIsActiveTr :: Salinity  is a dynamically active tracer
C     saltForcing    :: Flag which turns external forcing of salinity on/off
C     maskIniTemp    :: apply mask to initial Pot.Temp.
C     maskIniSalt    :: apply mask to initial salinity
C     checkIniTemp   :: check for points with identically zero initial Pot.Temp.
C     checkIniSalt   :: check for points with identically zero initial salinity
C- Pressure solver related parameters (PARM02)
C     useNSACGSolver :: Set to true to use "not self-adjoint" conjugate
C                       gradient solver that stores the iteration history
C                       for an iterative adjoint as accuate as possible
C     useSRCGSolver  :: Set to true to use conjugate gradient
C                       solver with single reduction (only one call of
C                       s/r mpi_allreduce), default is false
C- Time-stepping & free-surface params:
C     rigidLid            :: Set to true to use rigid lid
C     implicitFreeSurface :: Set to true to use implicit free surface
C     uniformLin_PhiSurf  :: Set to true to use a uniform Bo_surf in the
C                            linear relation Phi_surf = Bo_surf*eta
C     uniformFreeSurfLev  :: TRUE if free-surface level-index is uniform (=1)
C     exactConserv        :: Set to true to conserve exactly the total Volume
C     linFSConserveTr     :: Set to true to correct source/sink of tracer
C                            at the surface due to Linear Free Surface
C     useRealFreshWaterFlux :: if True (=Natural BCS), treats P+R-E flux
C                         as a real Fresh Water (=> changes the Sea Level)
C                         if F, converts P+R-E to salt flux (no SL effect)
C     storePhiHyd4Phys :: store hydrostatic potential for use in Physics/EOS
C                         this requires specific code for restart & exchange
C     quasiHydrostatic :: Using non-hydrostatic terms in hydrostatic algorithm
C     nonHydrostatic   :: Using non-hydrostatic algorithm
C     use3Dsolver      :: set to true to use 3-D pressure solver
C     implicitIntGravWave :: treat Internal Gravity Wave implicitly
C     staggerTimeStep   :: enable a Stagger time stepping U,V (& W) then T,S
C     applyExchUV_early :: Apply EXCH to U,V earlier, just before integr_continuity
C     doResetHFactors   :: Do reset thickness factors @ beginning of each time-step
C     implicitDiffusion :: Turns implicit vertical diffusion on
C     implicitViscosity :: Turns implicit vertical viscosity on
C     tempImplVertAdv   :: Turns on implicit vertical advection for Temperature
C     saltImplVertAdv   :: Turns on implicit vertical advection for Salinity
C     momImplVertAdv    :: Turns on implicit vertical advection for Momentum
C     multiDimAdvection :: Flag that enable multi-dimension advection
C     useMultiDimAdvec  :: True if multi-dim advection is used at least once
C     momDissip_In_AB   :: if False, put Dissipation tendency contribution
C                          out off Adams-Bashforth time stepping.
C     doAB_onGtGs       :: if the Adams-Bashforth time stepping is used, always
C                          apply AB on tracer tendencies (rather than on Tracer)
C- Other forcing params -
C     balanceEmPmR    :: substract global mean of EmPmR at every time step
C     balanceQnet     :: substract global mean of Qnet at every time step
C     balancePrintMean:: print substracted global means to STDOUT
C     doThetaClimRelax :: Set true if relaxation to temperature
C                        climatology is required.
C     doSaltClimRelax  :: Set true if relaxation to salinity
C                        climatology is required.
C     balanceThetaClimRelax :: substract global mean effect at every time step
C     balanceSaltClimRelax :: substract global mean effect at every time step
C     allowFreezing  :: Allows surface water to freeze and form ice
C     periodicExternalForcing :: Set true if forcing is time-dependant
C- I/O parameters -
C     globalFiles    :: Selects between "global" and "tiled" files.
C                       On some platforms with MPI, option globalFiles is either
C                       slow or does not work. Use useSingleCpuIO instead.
C     useSingleCpuIO :: moved to EEPARAMS.h
C     pickupStrictlyMatch :: check and stop if pickup-file do not stricly match
C     startFromPickupAB2 :: with AB-3 code, start from an AB-2 pickup
C     usePickupBeforeC54 :: start from old-pickup files, generated with code from
C                           before checkpoint-54a, Jul 06, 2004.
C     pickup_write_mdsio :: use mdsio to write pickups
C     pickup_read_mdsio  :: use mdsio to read  pickups
C     pickup_write_immed :: echo the pickup immediately (for conversion)
C     writePickupAtEnd   :: write pickup at the last timestep
C     timeave_mdsio      :: use mdsio for timeave output
C     snapshot_mdsio     :: use mdsio for "snapshot" (dumpfreq/diagfreq) output
C     monitor_stdio      :: use stdio for monitor output
C     dumpInitAndLast :: dumps model state to files at Initial (nIter0)
C                        & Last iteration, in addition multiple of dumpFreq iter.

      COMMON /PARM_L/
     & fluidIsAir, fluidIsWater,
     & usingPCoords, usingZCoords,
     & usingCartesianGrid, usingSphericalPolarGrid, rotateGrid,
     & usingCylindricalGrid, usingCurvilinearGrid, hasWetCSCorners,
     & deepAtmosphere, setInterFDr, setCenterDr, useMin4hFacEdges,
     & interViscAr_pCell, interDiffKr_pCell,
     & no_slip_sides, no_slip_bottom, bottomVisc_pCell, useSmag3D,
     & useFullLeith, useStrainTensionVisc, useAreaViscLength,
     & momViscosity, momAdvection, momForcing, momTidalForcing,
     & momPressureForcing, metricTerms, useNHMTerms,
     & useCoriolis, use3dCoriolis,
     & useCDscheme, vectorInvariantMomentum,
     & useJamartMomAdv, upwindVorticity, highOrderVorticity,
     & useAbsVorticity, upwindShear,
     & momStepping, calc_wVelocity, tempStepping, saltStepping,
     & addFrictionHeating, temp_stayPositive, salt_stayPositive,
     & tempAdvection, tempVertDiff4, tempIsActiveTr, tempForcing,
     & saltAdvection, saltVertDiff4, saltIsActiveTr, saltForcing,
     & maskIniTemp, maskIniSalt, checkIniTemp, checkIniSalt,
     & useNSACGSolver, useSRCGSolver,
     & rigidLid, implicitFreeSurface,
     & uniformLin_PhiSurf, uniformFreeSurfLev,
     & exactConserv, linFSConserveTr, useRealFreshWaterFlux,
     & storePhiHyd4Phys, quasiHydrostatic, nonHydrostatic,
     & use3Dsolver, implicitIntGravWave, staggerTimeStep,
     & applyExchUV_early, doResetHFactors,
     & implicitDiffusion, implicitViscosity,
     & tempImplVertAdv, saltImplVertAdv, momImplVertAdv,
     & multiDimAdvection, useMultiDimAdvec,
     & momDissip_In_AB, doAB_onGtGs,
     & balanceEmPmR, balanceQnet, balancePrintMean,
     & balanceThetaClimRelax, balanceSaltClimRelax,
     & doThetaClimRelax, doSaltClimRelax,
     & allowFreezing,
     & periodicExternalForcing,
     & globalFiles,
     & pickupStrictlyMatch, usePickupBeforeC54, startFromPickupAB2,
     & pickup_read_mdsio, pickup_write_mdsio, pickup_write_immed,
     & writePickupAtEnd,
     & timeave_mdsio, snapshot_mdsio, monitor_stdio,
     & outputTypesInclusive, dumpInitAndLast

      LOGICAL fluidIsAir
      LOGICAL fluidIsWater
      LOGICAL usingPCoords
      LOGICAL usingZCoords
      LOGICAL usingCartesianGrid
      LOGICAL usingSphericalPolarGrid, rotateGrid
      LOGICAL usingCylindricalGrid
      LOGICAL usingCurvilinearGrid, hasWetCSCorners
      LOGICAL deepAtmosphere
      LOGICAL setInterFDr
      LOGICAL setCenterDr
      LOGICAL useMin4hFacEdges
      LOGICAL interViscAr_pCell
      LOGICAL interDiffKr_pCell

      LOGICAL no_slip_sides
      LOGICAL no_slip_bottom
      LOGICAL bottomVisc_pCell
      LOGICAL useSmag3D
      LOGICAL useFullLeith
      LOGICAL useStrainTensionVisc
      LOGICAL useAreaViscLength
      LOGICAL momViscosity
      LOGICAL momAdvection
      LOGICAL momForcing
      LOGICAL momTidalForcing
      LOGICAL momPressureForcing
      LOGICAL metricTerms
      LOGICAL useNHMTerms

      LOGICAL useCoriolis
      LOGICAL use3dCoriolis
      LOGICAL useCDscheme
      LOGICAL vectorInvariantMomentum
      LOGICAL useJamartMomAdv
      LOGICAL upwindVorticity
      LOGICAL highOrderVorticity
      LOGICAL useAbsVorticity
      LOGICAL upwindShear
      LOGICAL momStepping
      LOGICAL calc_wVelocity
      LOGICAL tempStepping
      LOGICAL saltStepping
      LOGICAL addFrictionHeating
      LOGICAL temp_stayPositive
      LOGICAL salt_stayPositive
      LOGICAL tempAdvection
      LOGICAL tempVertDiff4
      LOGICAL tempIsActiveTr
      LOGICAL tempForcing
      LOGICAL saltAdvection
      LOGICAL saltVertDiff4
      LOGICAL saltIsActiveTr
      LOGICAL saltForcing
      LOGICAL maskIniTemp
      LOGICAL maskIniSalt
      LOGICAL checkIniTemp
      LOGICAL checkIniSalt
      LOGICAL useNSACGSolver
      LOGICAL useSRCGSolver
      LOGICAL rigidLid
      LOGICAL implicitFreeSurface
      LOGICAL uniformLin_PhiSurf
      LOGICAL uniformFreeSurfLev
      LOGICAL exactConserv
      LOGICAL linFSConserveTr
      LOGICAL useRealFreshWaterFlux
      LOGICAL storePhiHyd4Phys
      LOGICAL quasiHydrostatic
      LOGICAL nonHydrostatic
      LOGICAL use3Dsolver
      LOGICAL implicitIntGravWave
      LOGICAL staggerTimeStep
      LOGICAL applyExchUV_early
      LOGICAL doResetHFactors
      LOGICAL implicitDiffusion
      LOGICAL implicitViscosity
      LOGICAL tempImplVertAdv
      LOGICAL saltImplVertAdv
      LOGICAL momImplVertAdv
      LOGICAL multiDimAdvection
      LOGICAL useMultiDimAdvec
      LOGICAL momDissip_In_AB
      LOGICAL doAB_onGtGs
      LOGICAL balanceEmPmR
      LOGICAL balanceQnet
      LOGICAL balancePrintMean
      LOGICAL doThetaClimRelax
      LOGICAL doSaltClimRelax
      LOGICAL balanceThetaClimRelax
      LOGICAL balanceSaltClimRelax
      LOGICAL allowFreezing
      LOGICAL periodicExternalForcing
      LOGICAL globalFiles
      LOGICAL pickupStrictlyMatch
      LOGICAL usePickupBeforeC54
      LOGICAL startFromPickupAB2
      LOGICAL pickup_read_mdsio, pickup_write_mdsio
      LOGICAL pickup_write_immed, writePickupAtEnd
      LOGICAL timeave_mdsio, snapshot_mdsio, monitor_stdio
      LOGICAL outputTypesInclusive
      LOGICAL dumpInitAndLast

C--   COMMON /PARM_R/ "Real" valued parameters used by the model.
C     cg2dTargetResidual
C          :: Target residual for cg2d solver; no unit (RHS normalisation)
C     cg2dTargetResWunit
C          :: Target residual for cg2d solver; W unit (No RHS normalisation)
C     cg3dTargetResidual
C               :: Target residual for cg3d solver.
C     cg2dpcOffDFac :: Averaging weight for preconditioner off-diagonal.
C     Note. 20th May 1998
C           I made a weird discovery! In the model paper we argue
C           for the form of the preconditioner used here ( see
C           A Finite-volume, Incompressible Navier-Stokes Model
C           ...., Marshall et. al ). The algebra gives a simple
C           0.5 factor for the averaging of ac and aCw to get a
C           symmettric pre-conditioner. By using a factor of 0.51
C           i.e. scaling the off-diagonal terms in the
C           preconditioner down slightly I managed to get the
C           number of iterations for convergence in a test case to
C           drop form 192 -> 134! Need to investigate this further!
C           For now I have introduced a parameter cg2dpcOffDFac which
C           defaults to 0.51 but can be set at runtime.
C     delR      :: Vertical grid spacing ( units of r ).
C     delRc     :: Vertical grid spacing between cell centers (r unit).
C     delX      :: Separation between cell faces (m) or (deg), depending
C     delY         on input flags. Note: moved to header file SET_GRID.h
C     xgOrigin   :: Origin of the X-axis (Cartesian Grid) / Longitude of Western
C                :: most cell face (Lat-Lon grid) (Note: this is an "inert"
C                :: parameter but it makes geographical references simple.)
C     ygOrigin   :: Origin of the Y-axis (Cartesian Grid) / Latitude of Southern
C                :: most face (Lat-Lon grid).
C     rSphere    :: Radius of sphere for a spherical polar grid ( m ).
C     recip_rSphere :: Reciprocal radius of sphere ( m^-1 ).
C     radius_fromHorizGrid :: sphere Radius of input horiz. grid (Curvilinear Grid)
C     seaLev_Z   :: the reference height of sea-level (usually zero)
C     top_Pres   :: pressure (P-Coords) or reference pressure (Z-Coords) at the top
C     rSigmaBnd  :: vertical position (in r-unit) of r/sigma transition (Hybrid-Sigma)
C     gravity    :: Acceleration due to constant gravity ( m/s^2 )
C     recip_gravity :: Reciprocal gravity acceleration ( s^2/m )
C     gBaro      :: Accel. due to gravity used in barotropic equation ( m/s^2 )
C     gravFacC   :: gravity factor (vs surf. gravity) vert. profile at cell-Center
C     gravFacF   :: gravity factor (vs surf. gravity) vert. profile at cell-interF
C     rhoNil     :: Reference density for the linear equation of state
C     rhoConst   :: Vertically constant reference density (Boussinesq)
C     rho1Ref    :: reference vertical profile for density (anelastic)
C     rhoFacC    :: normalized (by rhoConst) reference density at cell-Center
C     rhoFacF    :: normalized (by rhoConst) reference density at cell-interFace
C     rhoConstFresh :: Constant reference density for fresh water (rain)
C     thetaConst :: Constant reference for potential temperature
C     tRef       :: reference vertical profile for potential temperature
C     sRef       :: reference vertical profile for salinity/specific humidity
C     surf_pRef  :: surface reference pressure ( Pa )
C     pRef4EOS   :: reference pressure used in EOS (case selectP_inEOS_Zc=1)
C     phiRef     :: reference potential (press/rho, geopot) profile (m^2/s^2)
C     dBdrRef    :: vertical gradient of reference buoyancy  [(m/s/r)^2]:
C                :: z-coord: = N^2_ref = Brunt-Vaissala frequency [s^-2]
C                :: p-coord: = -(d.alpha/dp)_ref          [(m^2.s/kg)^2]
C     rVel2wUnit :: units conversion factor (Non-Hydrostatic code),
C                :: from r-coordinate vertical velocity to vertical velocity [m/s].
C                :: z-coord: = 1 ; p-coord: wSpeed [m/s] = rVel [Pa/s] * rVel2wUnit
C     wUnit2rVel :: units conversion factor (Non-Hydrostatic code),
C                :: from vertical velocity [m/s] to r-coordinate vertical velocity.
C                :: z-coord: = 1 ; p-coord: rVel [Pa/s] = wSpeed [m/s] * wUnit2rVel
C     mass2rUnit :: units conversion factor (surface forcing),
C                :: from mass per unit area [kg/m2] to vertical r-coordinate unit.
C                :: z-coord: = 1/rhoConst ( [kg/m2] / rho = [m] ) ;
C                :: p-coord: = gravity    ( [kg/m2] *  g = [Pa] ) ;
C     rUnit2mass :: units conversion factor (surface forcing),
C                :: from vertical r-coordinate unit to mass per unit area [kg/m2].
C                :: z-coord: = rhoConst  ( [m] * rho = [kg/m2] ) ;
C                :: p-coord: = 1/gravity ( [Pa] /  g = [kg/m2] ) ;
C     sIceLoadFac:: factor to scale (and turn off) sIceLoad (sea-ice loading)
C                   default = 1
C     f0         :: Reference coriolis parameter ( 1/s )
C                   ( Southern edge f for beta plane )
C     beta       :: df/dy ( s^-1.m^-1 )
C     fPrime     :: Second Coriolis parameter ( 1/s ), related to Y-component
C                   of rotation (reference value = 2.Omega.Cos(Phi))
C     omega      :: Angular velocity ( rad/s )
C     rotationPeriod :: Rotation period (s) (= 2.pi/omega)
C     viscArNr   :: vertical profile of Eddy viscosity coeff.
C                   for vertical mixing of momentum ( units of r^2/s )
C     viscAh     :: Eddy viscosity coeff. for mixing of
C                   momentum laterally ( m^2/s )
C     viscAhW    :: Eddy viscosity coeff. for mixing of vertical
C                   momentum laterally, no effect for hydrostatic
C                   model, defaults to viscAhD if unset ( m^2/s )
C                   Not used if variable horiz. viscosity is used.
C     viscA4     :: Biharmonic viscosity coeff. for mixing of
C                   momentum laterally ( m^4/s )
C     viscA4W    :: Biharmonic viscosity coeff. for mixing of vertical
C                   momentum laterally, no effect for hydrostatic
C                   model, defaults to viscA4D if unset ( m^2/s )
C                   Not used if variable horiz. viscosity is used.
C     viscAhD    :: Eddy viscosity coeff. for mixing of momentum laterally
C                   (act on Divergence part) ( m^2/s )
C     viscAhZ    :: Eddy viscosity coeff. for mixing of momentum laterally
C                   (act on Vorticity  part) ( m^2/s )
C     viscA4D    :: Biharmonic viscosity coeff. for mixing of momentum laterally
C                   (act on Divergence part) ( m^4/s )
C     viscA4Z    :: Biharmonic viscosity coeff. for mixing of momentum laterally
C                   (act on Vorticity  part) ( m^4/s )
C     smag3D_coeff     :: Isotropic 3-D Smagorinsky viscosity coefficient (-)
C     smag3D_diffCoeff :: Isotropic 3-D Smagorinsky diffusivity coefficient (-)
C     viscC2leith  :: Leith non-dimensional viscosity factor (grad(vort))
C     viscC2leithD :: Modified Leith non-dimensional visc. factor (grad(div))
C     viscC2LeithQG:: QG Leith non-dimensional viscosity factor
C     viscC4leith  :: Leith non-dimensional viscosity factor (grad(vort))
C     viscC4leithD :: Modified Leith non-dimensional viscosity factor (grad(div))
C     viscC2smag   :: Smagorinsky non-dimensional viscosity factor (harmonic)
C     viscC4smag   :: Smagorinsky non-dimensional viscosity factor (biharmonic)
C     viscAhMax    :: Maximum eddy viscosity coeff. for mixing of
C                    momentum laterally ( m^2/s )
C     viscAhReMax  :: Maximum gridscale Reynolds number for eddy viscosity
C                     coeff. for mixing of momentum laterally (non-dim)
C     viscAhGrid   :: non-dimensional grid-size dependent viscosity
C     viscAhGridMax:: maximum and minimum harmonic viscosity coefficients ...
C     viscAhGridMin::  in terms of non-dimensional grid-size dependent visc.
C     viscA4Max    :: Maximum biharmonic viscosity coeff. for mixing of
C                     momentum laterally ( m^4/s )
C     viscA4ReMax  :: Maximum Gridscale Reynolds number for
C                     biharmonic viscosity coeff. momentum laterally (non-dim)
C     viscA4Grid   :: non-dimensional grid-size dependent bi-harmonic viscosity
C     viscA4GridMax:: maximum and minimum biharmonic viscosity coefficients ...
C     viscA4GridMin::  in terms of non-dimensional grid-size dependent viscosity
C     diffKhT   :: Laplacian diffusion coeff. for mixing of
C                 heat laterally ( m^2/s )
C     diffK4T   :: Biharmonic diffusion coeff. for mixing of
C                 heat laterally ( m^4/s )
C     diffKrNrT :: vertical profile of Laplacian diffusion coeff.
C                 for mixing of heat vertically ( units of r^2/s )
C     diffKr4T  :: vertical profile of Biharmonic diffusion coeff.
C                 for mixing of heat vertically ( units of r^4/s )
C     diffKhS  ::  Laplacian diffusion coeff. for mixing of
C                 salt laterally ( m^2/s )
C     diffK4S   :: Biharmonic diffusion coeff. for mixing of
C                 salt laterally ( m^4/s )
C     diffKrNrS :: vertical profile of Laplacian diffusion coeff.
C                 for mixing of salt vertically ( units of r^2/s ),
C     diffKr4S  :: vertical profile of Biharmonic diffusion coeff.
C                 for mixing of salt vertically ( units of r^4/s )
C     diffKrBL79surf :: T/S surface diffusivity (m^2/s) Bryan and Lewis, 1979
C     diffKrBL79deep :: T/S deep diffusivity (m^2/s) Bryan and Lewis, 1979
C     diffKrBL79scl  :: depth scale for arctan fn (m) Bryan and Lewis, 1979
C     diffKrBL79Ho   :: depth offset for arctan fn (m) Bryan and Lewis, 1979
C     BL79LatVary    :: polarwise of this latitude diffKrBL79 is applied with
C                       gradual transition to diffKrBLEQ towards Equator
C     diffKrBLEQsurf :: same as diffKrBL79surf but at Equator
C     diffKrBLEQdeep :: same as diffKrBL79deep but at Equator
C     diffKrBLEQscl  :: same as diffKrBL79scl but at Equator
C     diffKrBLEQHo   :: same as diffKrBL79Ho but at Equator
C     pCellMix_maxFac :: maximum enhanced mixing factor for thin partial-cell
C     pCellMix_delR   :: thickness criteria   for too thin partial-cell
C     pCellMix_viscAr :: vertical viscosity   for too thin partial-cell
C     pCellMix_diffKr :: vertical diffusivity for too thin partial-cell
C     deltaT    :: Default timestep ( s )
C     deltaTClock  :: Timestep used as model "clock". This determines the
C                    IO frequencies and is used in tagging output. It can
C                    be totally different to the dynamical time. Typically
C                    it will be the deep-water timestep for accelerated runs.
C                    Frequency of checkpointing and dumping of the model state
C                    are referenced to this clock. ( s )
C     deltaTMom    :: Timestep for momemtum equations ( s )
C     dTtracerLev  :: Timestep for tracer equations ( s ), function of level k
C     deltaTFreeSurf :: Timestep for free-surface equation ( s )
C     freeSurfFac  :: Parameter to turn implicit free surface term on or off
C                     freeSurFac = 1. uses implicit free surface
C                     freeSurFac = 0. uses rigid lid
C     abEps        :: Adams-Bashforth-2 stabilizing weight
C     alph_AB      :: Adams-Bashforth-3 primary factor
C     beta_AB      :: Adams-Bashforth-3 secondary factor
C     implicSurfPress :: parameter of the Crank-Nickelson time stepping :
C                     Implicit part of Surface Pressure Gradient ( 0-1 )
C     implicDiv2DFlow :: parameter of the Crank-Nickelson time stepping :
C                     Implicit part of barotropic flow Divergence ( 0-1 )
C     implicitNHPress :: parameter of the Crank-Nickelson time stepping :
C                     Implicit part of Non-Hydrostatic Pressure Gradient ( 0-1 )
C     hFacMin      :: Minimum fraction size of a cell (affects hFacC etc...)
C     hFacMinDz    :: Minimum dimensional size of a cell (affects hFacC etc..., m)
C     hFacMinDp    :: Minimum dimensional size of a cell (affects hFacC etc..., Pa)
C     hFacMinDr    :: Minimum dimensional size of a cell (-> hFacC etc..., r units)
C     hFacInf      :: Threshold (inf and sup) for fraction size of surface cell
C     hFacSup          that control vanishing and creating levels
C     tauCD         :: CD scheme coupling timescale ( s )
C     rCD           :: CD scheme normalised coupling parameter (= 1 - deltaT/tauCD)
C     epsAB_CD      :: Adams-Bashforth-2 stabilizing weight used in CD scheme
C     baseTime      :: model base time (time origin) = time @ iteration zero
C     startTime     :: Starting time for this integration ( s ).
C     endTime       :: Ending time for this integration ( s ).
C     chkPtFreq     :: Frequency of rolling check pointing ( s ).
C     pChkPtFreq    :: Frequency of permanent check pointing ( s ).
C     dumpFreq      :: Frequency with which model state is written to
C                      post-processing files ( s ).
C     diagFreq      :: Frequency with which model writes diagnostic output
C                      of intermediate quantities.
C     afFacMom      :: Advection of momentum term tracer parameter
C     vfFacMom      :: Momentum viscosity tracer parameter
C     pfFacMom      :: Momentum pressure forcing tracer parameter
C     cfFacMom      :: Coriolis term tracer parameter
C     foFacMom      :: Momentum forcing tracer parameter
C     mtFacMom      :: Metric terms tracer parameter
C     cosPower      :: Power of cosine of latitude to multiply viscosity
C     cAdjFreq      :: Frequency of convective adjustment
C
C     taveFreq      :: Frequency with which time-averaged model state
C                      is written to post-processing files ( s ).
C     tave_lastIter :: (for state variable only) fraction of the last time
C                      step (of each taveFreq period) put in the time average.
C                      (fraction for 1rst iter = 1 - tave_lastIter)
C     tauThetaClimRelax :: Relaxation to climatology time scale ( s ).
C     tauSaltClimRelax :: Relaxation to climatology time scale ( s ).
C     latBandClimRelax :: latitude band where Relaxation to Clim. is applied,
C                         i.e. where |yC| <= latBandClimRelax
C     externForcingPeriod :: Is the period of which forcing varies (eg. 1 month)
C     externForcingCycle :: Is the repeat time of the forcing (eg. 1 year)
C                          (note: externForcingCycle must be an integer
C                           number times externForcingPeriod)
C     convertFW2Salt :: salinity, used to convert Fresh-Water Flux to Salt Flux
C                       (use model surface (local) value if set to -1)
C     temp_EvPrRn :: temperature of Rain & Evap.
C     salt_EvPrRn :: salinity of Rain & Evap.
C     temp_addMass :: temperature of addMass array
C     salt_addMass :: salinity of addMass array
C        (notes: a) tracer content of Rain/Evap only used if both
C                     NonLin_FrSurf & useRealFreshWater are set.
C                b) use model surface (local) value if set to UNSET_RL)
C     hMixCriteria:: criteria for mixed-layer diagnostic
C     dRhoSmall   :: parameter for mixed-layer diagnostic
C     hMixSmooth  :: Smoothing parameter for mixed-layer diag (default=0=no smoothing)
C     ivdc_kappa  :: implicit vertical diffusivity for convection [m^2/s]
C     sideDragFactor     :: side-drag scaling factor (used only if no_slip_sides)
C                           (default=2: full drag ; =1: gives half-slip BC)
C     bottomDragLinear    :: Linear    bottom-drag coefficient (units of [r]/s)
C     bottomDragQuadratic :: Quadratic bottom-drag coefficient (units of [r]/m)
C               (if using zcoordinate, units becomes linear: m/s, quadratic: [-])
C     smoothAbsFuncRange :: 1/2 of interval around zero, for which FORTRAN ABS
C                           is to be replace by a smoother function
C                           (affects myabs, mymin, mymax)
C     nh_Am2        :: scales the non-hydrostatic terms and changes internal scales
C                      (i.e. allows convection at different Rayleigh numbers)
C     tCylIn        :: Temperature of the cylinder inner boundary
C     tCylOut       :: Temperature of the cylinder outer boundary
C     phiEuler      :: Euler angle, rotation about original z-axis
C     thetaEuler    :: Euler angle, rotation about new x-axis
C     psiEuler      :: Euler angle, rotation about new z-axis
      COMMON /PARM_R/ cg2dTargetResidual, cg2dTargetResWunit,
     & cg2dpcOffDFac, cg3dTargetResidual,
     & delR, delRc, xgOrigin, ygOrigin, rSphere, recip_rSphere,
     & radius_fromHorizGrid, seaLev_Z, top_Pres, rSigmaBnd,
     & deltaT, deltaTMom, dTtracerLev, deltaTFreeSurf, deltaTClock,
     & abEps, alph_AB, beta_AB,
     & f0, beta, fPrime, omega, rotationPeriod,
     & viscFacAdj, viscAh, viscAhW, smag3D_coeff, smag3D_diffCoeff,
     & viscAhMax, viscAhGrid, viscAhGridMax, viscAhGridMin,
     & viscC2leith, viscC2leithD, viscC2LeithQG,
     & viscC2smag, viscC4smag,
     & viscAhD, viscAhZ, viscA4D, viscA4Z,
     & viscA4, viscA4W, viscA4Max,
     & viscA4Grid, viscA4GridMax, viscA4GridMin,
     & viscAhReMax, viscA4ReMax,
     & viscC4leith, viscC4leithD, viscArNr,
     & diffKhT, diffK4T, diffKrNrT, diffKr4T,
     & diffKhS, diffK4S, diffKrNrS, diffKr4S,
     & diffKrBL79surf, diffKrBL79deep, diffKrBL79scl, diffKrBL79Ho,
     & BL79LatVary,
     & diffKrBLEQsurf, diffKrBLEQdeep, diffKrBLEQscl, diffKrBLEQHo,
     & pCellMix_maxFac, pCellMix_delR, pCellMix_viscAr, pCellMix_diffKr,
     & tauCD, rCD, epsAB_CD,
     & freeSurfFac, implicSurfPress, implicDiv2DFlow, implicitNHPress,
     & hFacMin, hFacMinDz, hFacInf, hFacSup,
     & gravity, recip_gravity, gBaro,
     & gravFacC, recip_gravFacC, gravFacF, recip_gravFacF,
     & rhoNil, rhoConst, recip_rhoConst, rho1Ref,
     & rhoFacC, recip_rhoFacC, rhoFacF, recip_rhoFacF, rhoConstFresh,
     & thetaConst, tRef, sRef, surf_pRef, pRef4EOS, phiRef, dBdrRef,
     & rVel2wUnit, wUnit2rVel, mass2rUnit, rUnit2mass,
     & baseTime, startTime, endTime,
     & chkPtFreq, pChkPtFreq, dumpFreq, adjDumpFreq,
     & diagFreq, taveFreq, tave_lastIter, monitorFreq, adjMonitorFreq,
     & afFacMom, vfFacMom, pfFacMom, cfFacMom, foFacMom, mtFacMom,
     & cosPower, cAdjFreq,
     & tauThetaClimRelax, tauSaltClimRelax, latBandClimRelax,
     & externForcingCycle, externForcingPeriod,
     & convertFW2Salt, temp_EvPrRn, salt_EvPrRn,
     & temp_addMass, salt_addMass, hFacMinDr, hFacMinDp,
     & ivdc_kappa, hMixCriteria, dRhoSmall, hMixSmooth,
     & sideDragFactor, bottomDragLinear, bottomDragQuadratic, nh_Am2,
     & smoothAbsFuncRange, sIceLoadFac,
     & tCylIn, tCylOut,
     & phiEuler, thetaEuler, psiEuler

      Real*8 cg2dTargetResidual
      Real*8 cg2dTargetResWunit
      Real*8 cg3dTargetResidual
      Real*8 cg2dpcOffDFac
      Real*8 delR(Nr)
      Real*8 delRc(Nr+1)
      Real*8 xgOrigin
      Real*8 ygOrigin
      Real*8 rSphere
      Real*8 recip_rSphere
      Real*8 radius_fromHorizGrid
      Real*8 seaLev_Z
      Real*8 top_Pres
      Real*8 rSigmaBnd
      Real*8 deltaT
      Real*8 deltaTClock
      Real*8 deltaTMom
      Real*8 dTtracerLev(Nr)
      Real*8 deltaTFreeSurf
      Real*8 abEps, alph_AB, beta_AB
      Real*8 f0
      Real*8 beta
      Real*8 fPrime
      Real*8 omega
      Real*8 rotationPeriod
      Real*8 freeSurfFac
      Real*8 implicSurfPress
      Real*8 implicDiv2DFlow
      Real*8 implicitNHPress
      Real*8 hFacMin
      Real*8 hFacMinDz
      Real*8 hFacMinDp
      Real*8 hFacMinDr
      Real*8 hFacInf
      Real*8 hFacSup
      Real*8 viscArNr(Nr)
      Real*8 viscFacAdj
      Real*8 viscAh
      Real*8 viscAhW
      Real*8 viscAhD
      Real*8 viscAhZ
      Real*8 smag3D_coeff, smag3D_diffCoeff
      Real*8 viscAhMax
      Real*8 viscAhReMax
      Real*8 viscAhGrid, viscAhGridMax, viscAhGridMin
      Real*8 viscC2leith
      Real*8 viscC2leithD
      Real*8 viscC2LeithQG
      Real*8 viscC2smag
      Real*8 viscA4
      Real*8 viscA4W
      Real*8 viscA4D
      Real*8 viscA4Z
      Real*8 viscA4Max
      Real*8 viscA4ReMax
      Real*8 viscA4Grid, viscA4GridMax, viscA4GridMin
      Real*8 viscC4leith
      Real*8 viscC4leithD
      Real*8 viscC4smag
      Real*8 diffKhT
      Real*8 diffK4T
      Real*8 diffKrNrT(Nr)
      Real*8 diffKr4T(Nr)
      Real*8 diffKhS
      Real*8 diffK4S
      Real*8 diffKrNrS(Nr)
      Real*8 diffKr4S(Nr)
      Real*8 diffKrBL79surf
      Real*8 diffKrBL79deep
      Real*8 diffKrBL79scl
      Real*8 diffKrBL79Ho
      Real*8 BL79LatVary
      Real*8 diffKrBLEQsurf
      Real*8 diffKrBLEQdeep
      Real*8 diffKrBLEQscl
      Real*8 diffKrBLEQHo
      Real*8 pCellMix_maxFac
      Real*8 pCellMix_delR
      Real*8 pCellMix_viscAr(Nr)
      Real*8 pCellMix_diffKr(Nr)
      Real*8 tauCD, rCD, epsAB_CD
      Real*8 gravity,       recip_gravity
      Real*8 gBaro
      Real*8 gravFacC(Nr),   recip_gravFacC(Nr)
      Real*8 gravFacF(Nr+1), recip_gravFacF(Nr+1)
      Real*8 rhoNil
      Real*8 rhoConst,      recip_rhoConst
      Real*8 rho1Ref(Nr)
      Real*8 rhoFacC(Nr),   recip_rhoFacC(Nr)
      Real*8 rhoFacF(Nr+1), recip_rhoFacF(Nr+1)
      Real*8 rhoConstFresh
      Real*8 thetaConst
      Real*8 tRef(Nr)
      Real*8 sRef(Nr)
      Real*8 surf_pRef, pRef4EOS(Nr)
      Real*8 phiRef(2*Nr+1)
      Real*8 dBdrRef(Nr)
      Real*8 rVel2wUnit(Nr+1), wUnit2rVel(Nr+1)
      Real*8 mass2rUnit, rUnit2mass
      Real*8 baseTime
      Real*8 startTime
      Real*8 endTime
      Real*8 chkPtFreq
      Real*8 pChkPtFreq
      Real*8 dumpFreq
      Real*8 adjDumpFreq
      Real*8 diagFreq
      Real*8 taveFreq
      Real*8 tave_lastIter
      Real*8 monitorFreq
      Real*8 adjMonitorFreq
      Real*8 afFacMom
      Real*8 vfFacMom
      Real*8 pfFacMom
      Real*8 cfFacMom
      Real*8 foFacMom
      Real*8 mtFacMom
      Real*8 cosPower
      Real*8 cAdjFreq
      Real*8 tauThetaClimRelax
      Real*8 tauSaltClimRelax
      Real*8 latBandClimRelax
      Real*8 externForcingCycle
      Real*8 externForcingPeriod
      Real*8 convertFW2Salt
      Real*8 temp_EvPrRn
      Real*8 salt_EvPrRn
      Real*8 temp_addMass
      Real*8 salt_addMass
      Real*8 ivdc_kappa
      Real*8 hMixCriteria
      Real*8 dRhoSmall
      Real*8 hMixSmooth
      Real*8 sideDragFactor
      Real*8 bottomDragLinear
      Real*8 bottomDragQuadratic
      Real*8 smoothAbsFuncRange
      Real*8 sIceLoadFac
      Real*8 nh_Am2
      Real*8 tCylIn, tCylOut
      Real*8 phiEuler, thetaEuler, psiEuler

C--   COMMON /PARM_A/ Thermodynamics constants ?
      COMMON /PARM_A/ HeatCapacity_Cp
      Real*8 HeatCapacity_Cp

C--   COMMON /PARM_ATM/ Atmospheric physical parameters (Ideal Gas EOS, ...)
C     celsius2K :: convert centigrade (Celsius) degree to Kelvin
C     atm_Po    :: standard reference pressure
C     atm_Cp    :: specific heat (Cp) of the (dry) air at constant pressure
C     atm_Rd    :: gas constant for dry air
C     atm_kappa :: kappa = R/Cp (R: constant of Ideal Gas EOS)
C     atm_Rq    :: water vapour specific volume anomaly relative to dry air
C                  (e.g. typical value = (29/18 -1) 10^-3 with q [g/kg])
C     integr_GeoPot :: option to select the way we integrate the geopotential
C                     (still a subject of discussions ...)
C     selectFindRoSurf :: select the way surf. ref. pressure (=Ro_surf) is
C             derived from the orography. Implemented: 0,1 (see INI_P_GROUND)
      COMMON /PARM_ATM/
     &            celsius2K,
     &            atm_Cp, atm_Rd, atm_kappa, atm_Rq, atm_Po,
     &            integr_GeoPot, selectFindRoSurf
      Real*8 celsius2K
      Real*8 atm_Po, atm_Cp, atm_Rd, atm_kappa, atm_Rq
      INTEGER integr_GeoPot, selectFindRoSurf

C----------------------------------------
C-- Logical flags for selecting packages
      LOGICAL useGAD
      LOGICAL useOBCS
      LOGICAL useSHAP_FILT
      LOGICAL useZONAL_FILT
      LOGICAL useOPPS
      LOGICAL usePP81
      LOGICAL useKL10
      LOGICAL useMY82
      LOGICAL useGGL90
      LOGICAL useKPP
      LOGICAL useGMRedi
      LOGICAL useDOWN_SLOPE
      LOGICAL useBBL
      LOGICAL useCAL
      LOGICAL useEXF
      LOGICAL useBulkForce
      LOGICAL useEBM
      LOGICAL useCheapAML
      LOGICAL useAUTODIFF
      LOGICAL useGrdchk
      LOGICAL useSMOOTH
      LOGICAL usePROFILES
      LOGICAL useECCO
      LOGICAL useCTRL
      LOGICAL useSBO
      LOGICAL useFLT
      LOGICAL usePTRACERS
      LOGICAL useGCHEM
      LOGICAL useRBCS
      LOGICAL useOffLine
      LOGICAL useMATRIX
      LOGICAL useFRAZIL
      LOGICAL useSEAICE
      LOGICAL useSALT_PLUME
      LOGICAL useShelfIce
      LOGICAL useStreamIce
      LOGICAL useICEFRONT
      LOGICAL useThSIce
      LOGICAL useLand
      LOGICAL useATM2d
      LOGICAL useAIM
      LOGICAL useAtm_Phys
      LOGICAL useFizhi
      LOGICAL useGridAlt
      LOGICAL useDiagnostics
      LOGICAL useREGRID
      LOGICAL useLayers
      LOGICAL useMNC
      LOGICAL useRunClock
      LOGICAL useEMBED_FILES
      LOGICAL useMYPACKAGE
      COMMON /PARM_PACKAGES/
     &        useGAD, useOBCS, useSHAP_FILT, useZONAL_FILT,
     &        useOPPS, usePP81, useKL10, useMY82, useGGL90, useKPP,
     &        useGMRedi, useBBL, useDOWN_SLOPE,
     &        useCAL, useEXF, useBulkForce, useEBM, useCheapAML,
     &        useGrdchk, useSMOOTH, usePROFILES, useECCO, useCTRL,
     &        useSBO, useFLT, useAUTODIFF,
     &        usePTRACERS, useGCHEM, useRBCS, useOffLine, useMATRIX,
     &        useFRAZIL, useSEAICE, useSALT_PLUME, useShelfIce,
     &        useStreamIce, useICEFRONT, useThSIce, useLand,
     &        useATM2D, useAIM, useAtm_Phys, useFizhi, useGridAlt,
     &        useDiagnostics, useREGRID, useLayers, useMNC,
     &        useRunClock, useEMBED_FILES,
     &        useMYPACKAGE

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
C     ==================================================================
C     HEADER exf_fields
C     ==================================================================
C
C     o Header file for the surface flux data.
C
C     started: Ralf.Giering@FastOpt.de 25-Mai-2000
C     changed: field swap in adj. mode; heimbach@mit.edu 10-Jan-2002
C     included runoff D. Stammer, Nov. 25, 2001
C     mods for pkg/seaice: menemenlis@jpl.nasa.gov 20-Dec-2002
C
C     ==================================================================
C     HEADER exf_fields
C     ==================================================================

C     Field definitions, units, and sign conventions:
C     ===============================================
C
C     ustress   :: Zonal surface wind stress in N/m^2
C                  > 0 for increase in uVel, which is west to
C                      east for cartesian and spherical polar grids
C                  Typical range: -0.5 < ustress < 0.5
C                  Input field
C
C     vstress   :: Meridional surface wind stress in N/m^2
C                  > 0 for increase in vVel, which is south to
C                      north for cartesian and spherical polar grids
C                  Typical range: -0.5 < vstress < 0.5
C                  Input field
C
C     hflux     :: Net upward surface heat flux including shortwave in W/m^2
C                  hflux = latent + sensible + lwflux + swflux
C                  > 0 for decrease in theta (ocean cooling)
C                  Typical range: -250 < hflux < 600
C                  Input field
C
C     sflux     :: Net upward freshwater flux in m/s
C                  sflux = evap - precip - runoff
C                  > 0 for increase in salt (ocean salinity)
C                  Typical range: -1e-7 < sflux < 1e-7
C                  Input field
C
C     swflux    :: Net upward shortwave radiation in W/m^2
C                  swflux = - ( swdown - ice and snow absorption - reflected )
C                  > 0 for decrease in theta (ocean cooling)
C                  Typical range: -350 < swflux < 0
C                  Input field
C
C     uwind     :: Surface (10-m) zonal wind velocity in m/s
C                  > 0 for increase in uVel, which is west to
C                      east for cartesian and spherical polar grids
C                  Typical range: -10 < uwind < 10
C                  Input or input/output field
C
C     vwind     :: Surface (10-m) meridional wind velocity in m/s
C                  > 0 for increase in vVel, which is south to
C                      north for cartesian and spherical polar grids
C                  Typical range: -10 < vwind < 10
C                  Input or input/output field
C
C     wspeed    :: Surface (10-m) wind speed in m/s
C                  >= 0 sqrt(u^2+v^2)
C                  Typical range: 0 < wspeed < 10
C                  Input or input/output field
C
C     atemp     :: Surface (2-m) air temperature in deg K
C                  Typical range: 200 < atemp < 300
C                  Input or input/output field
C
C     aqh       :: Surface (2m) specific humidity in kg/kg
C                  Typical range: 0 < aqh < 0.02
C                  Input or input/output field
C
C     hs        :: sensible heat flux into ocean in W/m^2
C                  > 0 for increase in theta (ocean warming)
C
C     hl        :: latent   heat flux into ocean in W/m^2
C                  > 0 for increase in theta (ocean warming)
C
C     lwflux    :: Net upward longwave radiation in W/m^2
C                  lwflux = - ( lwdown - ice and snow absorption - emitted )
C                  > 0 for decrease in theta (ocean cooling)
C                  Typical range: -20 < lwflux < 170
C                  Input field
C
C     evap      :: Evaporation in m/s
C                  > 0 for increase in salt (ocean salinity)
C                  Typical range: 0 < evap < 2.5e-7
C                  Input, input/output, or output field
C
C     precip    :: Total Precipitation (rain+snow) in m/s of liquid water
C                  > 0 for decrease in salt (ocean salinity)
C                  Typical range: 0 < precip < 5e-7
C                  Input or input/output field
C
C     snowprecip :: snow precipitation in m/s of equivalent liquid water
C                  > 0 for decrease in salt (ocean salinity)
C                  Typical range: 0 < precip < 5e-7
C                  Input or input/output field
C
C     runoff    :: River and glacier runoff in m/s
C                  > 0 for decrease in salt (ocean salinity)
C                  Typical range: 0 < runoff < ????
C                  Input or input/output field
C
C     runoftemp :: Temperature of runoff in deg C
C
C     saltflx   :: Net upward salt flux in (g/kg).kg/m^2/s = g/m^2/s
C                  > 0 for decrease in SSS.
C                  Typical origin: salty sea-ice formation / melting.
C
C     swdown    :: Downward shortwave radiation in W/m^2
C                  > 0 for increase in theta (ocean warming)
C                  Typical range: 0 < swdown < 450
C                  Input/output field
C
C     lwdown    :: Downward longwave radiation in W/m^2
C                  > 0 for increase in theta (ocean warming)
C                  Typical range: 50 < lwdown < 450
C                  Input/output field
C
C     apressure :: Atmospheric surface pressure field in Pa
C                  Typical range: 88000 < apressure < 108000
C                  Input field
C
C     tidePot   :: Tidal geopotential forcing in m^2/s^2
C                  Typical range: -10 < tidePot < +10
C                  Input field

C     NOTES:
C     ======
C
C     By default all surface forcing fields are defined at the center
C     of each grid (the rVel location in model/inc/GRID.h) unless
C     flags readStressOnAgrid or readStressOnCgrid are set.
C
C     Input and output units and sign conventions can be customized
C     using variables exf_inscal_* and exf_outscal_*, which are set
C     by exf_readparms.F
C
C     Output fields fu, fv, Qnet, Qsw, and EmPmR are
C     defined in FFIELDS.h
C
C     Arrays *0 and *1 below are used for temporal interpolation.
C

      COMMON /exf_stress_r/ ustress, vstress
      Real*8 ustress   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 vstress   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_ustress_r/ ustress0, ustress1
      Real*8 ustress0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 ustress1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_vstress_r/ vstress0, vstress1
      Real*8 vstress0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 vstress1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /exf_wspeed_r/ wspeed
      Real*8 wspeed   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_wspeed_r/ wspeed0, wspeed1
      Real*8 wspeed0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 wspeed1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /exf_atm_wind_r/ uwind, vwind
      Real*8 uwind     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 vwind     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_uwind_r/ uwind0, uwind1
      Real*8 uwind0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 uwind1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_vwind_r/ vwind0, vwind1
      Real*8 vwind0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 vwind1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /exf_netflux_r/ hflux, sflux
      Real*8 hflux     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 sflux     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_hflux_r/ hflux0, hflux1
      Real*8 hflux0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 hflux1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_sflux_r/ sflux0, sflux1
      Real*8 sflux0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 sflux1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /exf_atm_temp_r/ atemp, aqh, hs, hl, lwflux,
     &                        evap, precip, snowprecip
      Real*8 atemp     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 aqh       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 hs        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 hl        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 lwflux    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 evap      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 precip    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 snowprecip (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_atemp_r/ atemp0, atemp1
      Real*8 atemp0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 atemp1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_aqh_r/ aqh0, aqh1
      Real*8 aqh0      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 aqh1      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_lwflux_r/ lwflux0, lwflux1
      Real*8 lwflux0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 lwflux1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_precip_r/ precip0, precip1
      Real*8 precip0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 precip1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_snowprecip_r/ snowprecip0, snowprecip1
      Real*8 snowprecip0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 snowprecip1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C     wStress   :: wind-stress magnitude [Pa=N/m^2], @ grid-cell center
C     sh        :: wind-speed [m/s] (always larger than uMin)
      COMMON /exfl_wind_r/ wStress, cw, sw, sh
      Real*8 wStress   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 cw        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 sw        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 sh        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /exf_swflux_r/ swflux
      Real*8 swflux    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_swflux_r/ swflux0, swflux1
      Real*8 swflux0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 swflux1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /exf_rad_down_r/ swdown, lwdown
      Real*8 swdown    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 lwdown    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /exfl_rad_down_r/ swdown0, swdown1, lwdown0, lwdown1
      Real*8 swdown0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 swdown1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 lwdown0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 lwdown1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /exf_apressure_r/ apressure, apressure0, apressure1
      Real*8 apressure (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 apressure0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 apressure1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /exfl_runoff_r/ runoff, runoff0, runoff1
      Real*8 runoff    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 runoff0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 runoff1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


      COMMON /exfl_saltflx_r/ saltflx, saltflx0, saltflx1
      Real*8 saltflx   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 saltflx0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 saltflx1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)



      COMMON /exf_clim_sst_r/ climsst,
     &                        climsst0, climsst1
      Real*8 climsst       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 climsst0      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 climsst1      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /exf_clim_sss_r/ climsss,
     &                        climsss0, climsss1
      Real*8 climsss       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 climsss0      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 climsss1      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


C     !INPUT PARAMETERS:
C     siceFrac  :: seaice fraction [0-1]
C     diagGroup :: select which group of surface fluxes to process
C               :: =1 : wind-stress ; =2 : other surf. fluxes
C     myTime    :: current time in simulation [s]
C     myIter    :: current iteration number in simulation
C     myThid    :: my Thread Id number
      Real*8     siceFrac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER diagGroup
      Real*8     myTime
      INTEGER myIter
      INTEGER myThid

C     !FUNCTIONS:
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON

C     !LOCAL VARIABLES:
C     tmpFld    :: open-ocean weighted surface flux
C     bi,bj     :: tile indices
C     i,j       :: grid point indices
      Real*8     tmpFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER bi, bj
      INTEGER i, j
CEOP

      IF ( useDiagnostics .AND. diags_opOceWeighted ) THEN

C--   First group of surface-flux: Momentum fluxes
        IF ( diagGroup.EQ.-1 ) THEN
          CALL DIAGNOSTICS_FILL( ustress,'EXFtaux ', 0,1,0,1,1,myThid )
          CALL DIAGNOSTICS_FILL( vstress,'EXFtauy ', 0,1,0,1,1,myThid )
        ELSEIF ( diagGroup.EQ.1 ) THEN

         IF ( DIAGNOSTICS_IS_ON('EXFtaux ',myThid) ) THEN
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
            IF ( stressIsOnCgrid ) THEN
             DO j = 1,sNy
              DO i = 1,sNx+1
               tmpFld(i,j,bi,bj) = ustress(i,j,bi,bj)
     &            *( oneRL - halfRL*
     &               ( siceFrac(i-1,j,bi,bj) + siceFrac(i,j,bi,bj) ) )
              ENDDO
             ENDDO
            ELSE
             DO j = 1-OLy,sNy+OLy
              DO i = 1-OLx,sNx+OLx
               tmpFld(i,j,bi,bj) = ustress(i,j,bi,bj)
     &                           *( oneRL - siceFrac(i,j,bi,bj) )
              ENDDO
             ENDDO
            ENDIF
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmpFld, 'EXFtaux ', 0,1,0,1,1,myThid )
         ENDIF

         IF ( DIAGNOSTICS_IS_ON('EXFtauy ',myThid) ) THEN
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
            IF ( stressIsOnCgrid ) THEN
             DO j = 1,sNy+1
              DO i = 1,sNx
               tmpFld(i,j,bi,bj) = vstress(i,j,bi,bj)
     &            *( oneRL - halfRL*
     &               ( siceFrac(i,j-1,bi,bj) + siceFrac(i,j,bi,bj) ) )
              ENDDO
             ENDDO
            ELSE
             DO j = 1-OLy,sNy+OLy
              DO i = 1-OLx,sNx+OLx
               tmpFld(i,j,bi,bj) = vstress(i,j,bi,bj)
     &                           *( oneRL - siceFrac(i,j,bi,bj) )
              ENDDO
             ENDDO
            ENDIF
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmpFld, 'EXFtauy ', 0,1,0,1,1,myThid )
         ENDIF

C--   end first group
        ENDIF

C--   Second group of surface-flux: Buoyancy fluxes
        IF ( diagGroup.EQ.-2 ) THEN
          CALL DIAGNOSTICS_FILL( hflux,  'EXFqnet ', 0,1,0,1,1,myThid )
          CALL DIAGNOSTICS_FILL( sflux,  'EXFempmr', 0,1,0,1,1,myThid )
          CALL DIAGNOSTICS_FILL( swflux, 'EXFswnet', 0,1,0,1,1,myThid )
        ELSEIF ( diagGroup.EQ.2 ) THEN

         IF ( DIAGNOSTICS_IS_ON('EXFqnet ',myThid) ) THEN
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
             DO j = 1-OLy,sNy+OLy
              DO i = 1-OLx,sNx+OLx
               tmpFld(i,j,bi,bj) = hflux(i,j,bi,bj)
     &                           *( oneRL - siceFrac(i,j,bi,bj) )
              ENDDO
             ENDDO
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmpFld, 'EXFqnet ', 0,1,0,1,1,myThid )
         ENDIF

         IF ( DIAGNOSTICS_IS_ON('EXFempmr',myThid) ) THEN
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
             DO j = 1-OLy,sNy+OLy
              DO i = 1-OLx,sNx+OLx
               tmpFld(i,j,bi,bj) = sflux(i,j,bi,bj)
     &                           *( oneRL - siceFrac(i,j,bi,bj) )
              ENDDO
             ENDDO
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmpFld, 'EXFempmr', 0,1,0,1,1,myThid )
         ENDIF

         IF ( DIAGNOSTICS_IS_ON('EXFswnet',myThid) ) THEN
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
             DO j = 1-OLy,sNy+OLy
              DO i = 1-OLx,sNx+OLx
               tmpFld(i,j,bi,bj) = swflux(i,j,bi,bj)
     &                           *( oneRL - siceFrac(i,j,bi,bj) )
              ENDDO
             ENDDO
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmpFld, 'EXFswnet', 0,1,0,1,1,myThid )
         ENDIF

C--   end second group
        ENDIF

C--   Third group: Individual components of surface buoyancy fluxes ;
C     processed the same way and at same time as buoyancy fluxes (second group)
        IF ( diagGroup.EQ.-2 ) THEN
          CALL DIAGNOSTICS_FILL( hs    , 'EXFhs   ', 0,1,0,1,1,myThid )
          CALL DIAGNOSTICS_FILL( hl    , 'EXFhl   ', 0,1,0,1,1,myThid )
          CALL DIAGNOSTICS_FILL( lwflux, 'EXFlwnet', 0,1,0,1,1,myThid )
          CALL DIAGNOSTICS_FILL( evap,   'EXFevap ', 0,1,0,1,1,myThid )
        ELSEIF ( diagGroup.EQ.2 ) THEN

         IF ( DIAGNOSTICS_IS_ON('EXFhs   ',myThid) ) THEN
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
             DO j = 1-OLy,sNy+OLy
              DO i = 1-OLx,sNx+OLx
               tmpFld(i,j,bi,bj) = hs(i,j,bi,bj)
     &                           *( oneRL - siceFrac(i,j,bi,bj) )
              ENDDO
             ENDDO
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmpFld, 'EXFhs   ', 0,1,0,1,1,myThid )
         ENDIF

         IF ( DIAGNOSTICS_IS_ON('EXFhl   ',myThid) ) THEN
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
             DO j = 1-OLy,sNy+OLy
              DO i = 1-OLx,sNx+OLx
               tmpFld(i,j,bi,bj) = hl(i,j,bi,bj)
     &                           *( oneRL - siceFrac(i,j,bi,bj) )
              ENDDO
             ENDDO
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmpFld, 'EXFhl   ', 0,1,0,1,1,myThid )
         ENDIF

         IF ( DIAGNOSTICS_IS_ON('EXFlwnet',myThid) ) THEN
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
             DO j = 1-OLy,sNy+OLy
              DO i = 1-OLx,sNx+OLx
               tmpFld(i,j,bi,bj) = lwflux(i,j,bi,bj)
     &                           *( oneRL - siceFrac(i,j,bi,bj) )
              ENDDO
             ENDDO
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmpFld, 'EXFlwnet', 0,1,0,1,1,myThid )
         ENDIF

         IF ( DIAGNOSTICS_IS_ON('EXFevap ',myThid) ) THEN
          DO bj=myByLo(myThid),myByHi(myThid)
           DO bi=myBxLo(myThid),myBxHi(myThid)
             DO j = 1-OLy,sNy+OLy
              DO i = 1-OLx,sNx+OLx
               tmpFld(i,j,bi,bj) = evap(i,j,bi,bj)
     &                           *( oneRL - siceFrac(i,j,bi,bj) )
              ENDDO
             ENDDO
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmpFld, 'EXFevap ', 0,1,0,1,1,myThid )
         ENDIF

C--   end third group
        ENDIF

      ENDIF

      RETURN
      END
