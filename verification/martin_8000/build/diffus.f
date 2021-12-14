












C     *==========================================================*
C     | SEAICE_OPTIONS.h
C     | o CPP options file for sea ice package.
C     *==========================================================*
C     | Use this file for selecting options within the sea ice
C     | package.
C     *==========================================================*








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

C--   Write "text-plots" of certain fields in STDOUT for debugging.

C--   Allow sea-ice dynamic code.
C     This option is provided to allow use of TAMC
C     on the thermodynamics component of the code only.
C     Sea-ice dynamics can also be turned off at runtime
C     using variable SEAICEuseDYNAMICS.

C--   By default, the sea-ice package uses its own integrated bulk
C     formulae to compute fluxes (fu, fv, EmPmR, Qnet, and Qsw) over
C     open-ocean.  When this flag is set, these variables are computed
C     in a separate external package, for example, pkg/exf, and then
C     modified for sea-ice effects by pkg/seaice.

C--   This CPP flag has been retired.  The number of ice categories
C     used to solve for seaice flux is now specified by run-time
C     parameter SEAICE_multDim.
C     Note: be aware of pickup_seaice.* compatibility issues when
C     restarting a simulation with a different number of categories.
c#define SEAICE_MULTICATEGORY

C--   run with sea Ice Thickness Distribution (ITD);
C     set number of categories (nITD) in SEAICE_SIZE.h

C--   Since the missing sublimation term is now included
C     this flag is needed for backward compatibility

C--   Suspected missing term in coupled ocn-ice heat budget (to be confirmed)

C--   Default is constant seaice salinity (SEAICE_salt0); Define the following
C     flag to consider (space & time) variable salinity: advected and forming
C     seaice with a fraction (=SEAICE_saltFrac) of freezing seawater salinity.
C- Note: SItracer also offers an alternative way to handle variable salinity.

C--   Tracers of ice and/or ice cover.

C--   Enable grease ice parameterization
C     The grease ice parameterization delays formation of solid
C     sea ice from frazil ice by a time constant and provides a
C     dynamic calculation of the initial solid sea ice thickness
C     HO as a function of winds, currents and available grease ice
C     volume. Grease ice does not significantly reduce heat loss
C     from the ocean in winter and area covered by grease is thus
C     handled like open water.
C     (For details see Smedsrud and Martin, 2014, Ann.Glac.)
C     Set SItrName(1) = 'grease' in namelist SEAICE_PARM03 in data.seaice
C     then output SItr01 is SItrNameLong(1) = 'grease ice volume fraction',
C     with SItrUnit(1) = '[0-1]', which needs to be multiplied by SIheff
C     to yield grease ice volume. Additionally, the actual grease ice
C     layer thickness (diagnostic SIgrsLT) can be saved.
C--   grease ice uses SItracer:

C--   Historically, the seaice model was discretized on a B-Grid. This
C     discretization should still work but it is not longer actively tested
C     and supported. The following flag should always be set in order to use
C     the operational C-grid discretization.

C--   Only for the C-grid version it is possible to
C     enable advection of sea ice momentum
C     enable JFNK code by defining the following flag
C     enable Krylov code by defining the following flag
C     enable this flag to reproduce old verification results for JFNK
C     enable LSR to use global (multi-tile) tri-diagonal solver
C     enable EVP code by defining the following flag
C--   When set use SEAICE_zetaMin and SEAICE_evpDampC to limit viscosities
C     from below and above in seaice_evp: not necessary, and not recommended
C     smooth regularization (without max-function) of delta for
C     better differentiability
C     regularize zeta to zmax with a smooth tanh-function instead
C     of a min(zeta,zmax). This improves convergence of iterative
C     solvers (Lemieux and Tremblay 2009, JGR). No effect on EVP
C     allow the truncated ellipse rheology (runtime flag SEAICEuseTEM)
C     allow using a damage parameter that records violations of a Mohr
C     Coulomb criterion for stress
C     allow Maxwell Elasto-Brittle rheology (runtime flat SEAICEuseMEB);
C     the MEB code re-uses a lot of the VP code and can be used together
C     with all implicit solvers, although the parameter tuning requires
C     some care and not all combinations give satisfactory results; you need
C     an accurate solution of the momentum equations to reduce the amount
C     of error accumulation
C     MEB always needs the damage field
C     Use LSR vector code; not useful on non-vector machines, because it
C     slows down convergence considerably, but the extra iterations are
C     more than made up by the much faster code on vector machines. For
C     the only regularly test vector machine these flags a specified
C     in the build options file SUPER-UX_SX-8_sxf90_awi, so that we comment
C     them out here.
C     Use zebra-method (alternate lines) for line-successive-relaxation
C     This modification improves the convergence of the vector code
C     dramatically, so that is may actually be useful in general, but
C     that needs to be tested. Can be used without vectorization options.
C     Use parameterisation of grounding ice for a better representation
C     of fastice in shallow seas

C--   When set limit the Ice-Loading to mass of 1/5 of Surface ocean grid-box
C--   When set use SEAICE_clipVelocties = .true., to clip U/VICE at 40cm/s,
C     not recommended
C--   When set cap the sublimation latent heat flux in solve4temp according
C     to the available amount of ice+snow. Otherwise this term is treated
C     like all of the others -- residuals heat and fw stocks are passed to
C     the ocean at the end of seaice_growth in a conservative manner.
C     SEAICE_CAP_SUBLIM is not needed as of now, but kept just in case.

C--   Enable the adjointable sea-ice thermodynamic model
C     uses seaice_growth_adx.F and seaice_solve4temp_adx.F

C--   Enable free drift code

C--   pkg/seaice cost functions compile flags
c       >>> Sea-ice volume (requires pkg/cost)
c       >>> Sea-ice misfit to obs (requires pkg/cost and ecco)


CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***

CBOP
C     !ROUTINE: DIFFUS
C     !INTERFACE:
      SUBROUTINE DIFFUS(
     U                   fld,
     I                   DIFFA, iceMsk, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R DIFFUS
C     | o Calculate laplacian of input field
C     |   and return the result in the same array
C     *==========================================================*
C     \ev

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
C
CBOP
C    !ROUTINE: GRID.h
C    !INTERFACE:
C    include GRID.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | GRID.h
C     | o Header file defining model grid.
C     *==========================================================*
C     | Model grid is defined for each process by reference to
C     | the arrays set here.
C     | Notes
C     | =====
C     | The standard MITgcm convention of westmost, southern most
C     | and upper most having the (1,1,1) index is used here.
C     | i.e.
C     |----------------------------------------------------------
C     | (1)  Plan view schematic of model grid (top layer i.e. )
C     |      ================================= ( ocean surface )
C     |                                        ( or top of     )
C     |                                        ( atmosphere    )
C     |      This diagram shows the location of the model
C     |      prognostic variables on the model grid. The "T"
C     |      location is used for all tracers. The figure also
C     |      shows the southern most, western most indexing
C     |      convention that is used for all model variables.
C     |
C     |
C     |             V(i=1,                     V(i=Nx,
C     |               j=Ny+1,                    j=Ny+1,
C     |               k=1)                       k=1)
C     |                /|\                       /|\  "PWX"
C     |       |---------|------------------etc..  |---- *---
C     |       |                     |                   *  |
C     |"PWY"*******************************etc..  **********"PWY"
C     |       |                     |                   *  |
C     |       |                     |                   *  |
C     |       |                     |                   *  |
C     |U(i=1, ==>       x           |             x     *==>U
C     |  j=Ny,|      T(i=1,         |          T(i=Nx,  *(i=Nx+1,
C     |  k=1) |        j=Ny,        |            j=Ny,  *  |j=Ny,
C     |       |        k=1)         |            k=1)   *  |k=1)
C     |
C     |       .                     .                      .
C     |       .                     .                      .
C     |       .                     .                      .
C     |       e                     e                   *  e
C     |       t                     t                   *  t
C     |       c                     c                   *  c
C     |       |                     |                   *  |
C     |       |                     |                   *  |
C     |U(i=1, ==>       x           |             x     *  |
C     |  j=2, |      T(i=1,         |          T(i=Nx,  *  |
C     |  k=1) |        j=2,         |            j=2,   *  |
C     |       |        k=1)         |            k=1)   *  |
C     |       |                     |                   *  |
C     |       |        /|\          |            /|\    *  |
C     |      -----------|------------------etc..  |-----*---
C     |       |       V(i=1,        |           V(i=Nx, *  |
C     |       |         j=2,        |             j=2,  *  |
C     |       |         k=1)        |             k=1)  *  |
C     |       |                     |                   *  |
C     |U(i=1, ==>       x         ==>U(i=2,       x     *==>U
C     |  j=1, |      T(i=1,         |  j=1,    T(i=Nx,  *(i=Nx+1,
C     |  k=1) |        j=1,         |  k=1)      j=1,   *  |j=1,
C     |       |        k=1)         |            k=1)   *  |k=1)
C     |       |                     |                   *  |
C     |       |        /|\          |            /|\    *  |
C     |"SB"++>|---------|------------------etc..  |-----*---
C     |      /+\      V(i=1,                    V(i=Nx, *
C     |       +         j=1,                      j=1,  *
C     |       +         k=1)                      k=1)  *
C     |     "WB"                                      "PWX"
C     |
C     |   N, y increasing northwards
C     |  /|\ j increasing northwards
C     |   |
C     |   |
C     |   ======>E, x increasing eastwards
C     |             i increasing eastwards
C     |
C     |    i: East-west index
C     |    j: North-south index
C     |    k: up-down index
C     |    U: x-velocity (m/s)
C     |    V: y-velocity (m/s)
C     |    T: potential temperature (oC)
C     | "SB": Southern boundary
C     | "WB": Western boundary
C     |"PWX": Periodic wrap around in X.
C     |"PWY": Periodic wrap around in Y.
C     |----------------------------------------------------------
C     | (2) South elevation schematic of model grid
C     |     =======================================
C     |     This diagram shows the location of the model
C     |     prognostic variables on the model grid. The "T"
C     |     location is used for all tracers. The figure also
C     |     shows the upper most, western most indexing
C     |     convention that is used for all model variables.
C     |
C     |      "WB"
C     |       +
C     |       +
C     |      \+/       /|\                       /|\       .
C     |"UB"++>|-------- | -----------------etc..  | ----*---
C     |       |    rVel(i=1,        |        rVel(i=Nx, *  |
C     |       |         j=1,        |             j=1,  *  |
C     |       |         k=1)        |             k=1)  *  |
C     |       |                     |                   *  |
C     |U(i=1, ==>       x         ==>U(i=2,       x     *==>U
C     |  j=1, |      T(i=1,         |  j=1,    T(i=Nx,  *(i=Nx+1,
C     |  k=1) |        j=1,         |  k=1)      j=1,   *  |j=1,
C     |       |        k=1)         |            k=1)   *  |k=1)
C     |       |                     |                   *  |
C     |       |        /|\          |            /|\    *  |
C     |       |-------- | -----------------etc..  | ----*---
C     |       |    rVel(i=1,        |        rVel(i=Nx, *  |
C     |       |         j=1,        |             j=1,  *  |
C     |       |         k=2)        |             k=2)  *  |
C     |
C     |       .                     .                      .
C     |       .                     .                      .
C     |       .                     .                      .
C     |       e                     e                   *  e
C     |       t                     t                   *  t
C     |       c                     c                   *  c
C     |       |                     |                   *  |
C     |       |                     |                   *  |
C     |       |                     |                   *  |
C     |       |                     |                   *  |
C     |       |        /|\          |            /|\    *  |
C     |       |-------- | -----------------etc..  | ----*---
C     |       |    rVel(i=1,        |        rVel(i=Nx, *  |
C     |       |         j=1,        |             j=1,  *  |
C     |       |         k=Nr)       |             k=Nr) *  |
C     |U(i=1, ==>       x         ==>U(i=2,       x     *==>U
C     |  j=1, |      T(i=1,         |  j=1,    T(i=Nx,  *(i=Nx+1,
C     |  k=Nr)|        j=1,         |  k=Nr)     j=1,   *  |j=1,
C     |       |        k=Nr)        |            k=Nr)  *  |k=Nr)
C     |       |                     |                   *  |
C     |"LB"++>==============================================
C     |                                               "PWX"
C     |
C     | Up   increasing upwards.
C     |/|\                                                       .
C     | |
C     | |
C     | =====> E  i increasing eastwards
C     | |         x increasing eastwards
C     | |
C     |\|/
C     | Down,k increasing downwards.
C     |
C     | Note: r => height (m) => r increases upwards
C     |       r => pressure (Pa) => r increases downwards
C     |
C     |
C     |    i: East-west index
C     |    j: North-south index
C     |    k: up-down index
C     |    U: x-velocity (m/s)
C     | rVel: z-velocity ( units of r )
C     |       The vertical velocity variable rVel is in units of
C     |       "r" the vertical coordinate. r in m will give
C     |       rVel m/s. r in Pa will give rVel Pa/s.
C     |    T: potential temperature (oC)
C     | "UB": Upper boundary.
C     | "LB": Lower boundary (always solid - therefore om|w == 0)
C     | "WB": Western boundary
C     |"PWX": Periodic wrap around in X.
C     |----------------------------------------------------------
C     | (3) Views showing nomenclature and indexing
C     |     for grid descriptor variables.
C     |
C     |      Fig 3a. shows the orientation, indexing and
C     |      notation for the grid spacing terms used internally
C     |      for the evaluation of gradient and averaging terms.
C     |      These varaibles are set based on the model input
C     |      parameters which define the model grid in terms of
C     |      spacing in X, Y and Z.
C     |
C     |      Fig 3b. shows the orientation, indexing and
C     |      notation for the variables that are used to define
C     |      the model grid. These varaibles are set directly
C     |      from the model input.
C     |
C     | Figure 3a
C     | =========
C     |       |------------------------------------
C     |       |                       |
C     |"PWY"********************************* etc...
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |
C     |       .                       .
C     |       .                       .
C     |       .                       .
C     |       e                       e
C     |       t                       t
C     |       c                       c
C     |       |-----------v-----------|-----------v----------|-
C     |       |                       |                      |
C     |       |                       |                      |
C     |       |                       |                      |
C     |       |                       |                      |
C     |       |                       |                      |
C     |       u<--dxF(i=1,j=2,k=1)--->u           t          |
C     |       |/|\       /|\          |                      |
C     |       | |         |           |                      |
C     |       | |         |           |                      |
C     |       | |         |           |                      |
C     |       |dyU(i=1,  dyC(i=1,     |                      |
C     | ---  ---|--j=2,---|--j=2,-----------------v----------|-
C     | /|\   | |  k=1)   |  k=1)     |          /|\         |
C     |  |    | |         |           |          dyF(i=2,    |
C     |  |    | |         |           |           |  j=1,    |
C     |dyG(   |\|/       \|/          |           |  k=1)    |
C     |   i=1,u---        t<---dxC(i=2,j=1,k=1)-->t          |
C     |   j=1,|                       |           |          |
C     |   k=1)|                       |           |          |
C     |  |    |                       |           |          |
C     |  |    |                       |           |          |
C     | \|/   |           |<---dxV(i=2,j=1,k=1)--\|/         |
C     |"SB"++>|___________v___________|___________v__________|_
C     |       <--dxG(i=1,j=1,k=1)----->
C     |      /+\                                              .
C     |       +
C     |       +
C     |     "WB"
C     |
C     |   N, y increasing northwards
C     |  /|\ j increasing northwards
C     |   |
C     |   |
C     |   ======>E, x increasing eastwards
C     |             i increasing eastwards
C     |
C     |    i: East-west index
C     |    j: North-south index
C     |    k: up-down index
C     |    u: x-velocity point
C     |    V: y-velocity point
C     |    t: tracer point
C     | "SB": Southern boundary
C     | "WB": Western boundary
C     |"PWX": Periodic wrap around in X.
C     |"PWY": Periodic wrap around in Y.
C     |
C     | Figure 3b
C     | =========
C     |
C     |       .                       .
C     |       .                       .
C     |       .                       .
C     |       e                       e
C     |       t                       t
C     |       c                       c
C     |       |-----------v-----------|-----------v--etc...
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       u<--delX(i=1)---------->u           t
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       |                       |
C     |       |-----------v-----------------------v--etc...
C     |       |          /|\          |
C     |       |           |           |
C     |       |           |           |
C     |       |           |           |
C     |       u        delY(j=1)      |           t
C     |       |           |           |
C     |       |           |           |
C     |       |           |           |
C     |       |           |           |
C     |       |          \|/          |
C     |"SB"++>|___________v___________|___________v__etc...
C     |      /+\                                                 .
C     |       +
C     |       +
C     |     "WB"
C     |
C     *==========================================================*
C     \ev
CEOP

C     Macros that override/modify standard definitions
C
CBOP
C    !ROUTINE: GRID_MACROS.h
C    !INTERFACE:
C    include GRID_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | GRID_MACROS.h
C     *==========================================================*
C     | These macros are used to substitute definitions for
C     | GRID.h variables for particular configurations.
C     | In setting these variables the following convention
C     | applies.
C     | undef  phi_CONST   - Indicates the variable phi is fixed
C     |                      in X, Y and Z.
C     | undef  phi_FX      - Indicates the variable phi only
C     |                      varies in X (i.e.not in X or Z).
C     | undef  phi_FY      - Indicates the variable phi only
C     |                      varies in Y (i.e.not in X or Z).
C     | undef  phi_FXY     - Indicates the variable phi only
C     |                      varies in X and Y ( i.e. not Z).
C     *==========================================================*
C     \ev
CEOP

C
CBOP
C    !ROUTINE: DXC_MACROS.h
C    !INTERFACE:
C    include DXC_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | DXC_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: DXF_MACROS.h
C    !INTERFACE:
C    include DXF_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | DXF_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: DXG_MACROS.h
C    !INTERFACE:
C    include DXG_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | DXG_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: DXV_MACROS.h
C    !INTERFACE:
C    include DXV_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | DXV_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: DYC_MACROS.h
C    !INTERFACE:
C    include DYC_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | DYC_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: DYF_MACROS.h
C    !INTERFACE:
C    include DYF_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | DYF_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: DYG_MACROS.h
C    !INTERFACE:
C    include DYG_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | DYG_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: DYU_MACROS.h
C    !INTERFACE:
C    include DYU_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | DYU_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: HFACC_MACROS.h
C    !INTERFACE:
C    include HFACC_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | HFACC_MACROS.h                                            
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP







C
CBOP
C    !ROUTINE: HFACS_MACROS.h
C    !INTERFACE:
C    include HFACS_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | HFACS_MACROS.h                                            
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP







C
CBOP
C    !ROUTINE: HFACW_MACROS.h
C    !INTERFACE:
C    include HFACW_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | HFACW_MACROS.h                                            
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP







C
CBOP
C    !ROUTINE: RECIP_DXC_MACROS.h
C    !INTERFACE:
C    include RECIP_DXC_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_DXC_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: RECIP_DXF_MACROS.h
C    !INTERFACE:
C    include RECIP_DXF_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_DXF_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: RECIP_DXG_MACROS.h
C    !INTERFACE:
C    include RECIP_DXG_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_DXG_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: RECIP_DXV_MACROS.h
C    !INTERFACE:
C    include RECIP_DXV_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_DXV_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: RECIP_DYC_MACROS.h
C    !INTERFACE:
C    include RECIP_DYC_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_DYC_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: RECIP_DYF_MACROS.h
C    !INTERFACE:
C    include RECIP_DYF_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_DYF_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: RECIP_DYG_MACROS.h
C    !INTERFACE:
C    include RECIP_DYG_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_DYG_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: RECIP_DYU_MACROS.h
C    !INTERFACE:
C    include RECIP_DYU_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_DYU_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: RECIP_HFACC_MACROS.h
C    !INTERFACE:
C    include RECIP_HFACC_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_HFACC_MACROS.h                                      
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP







C
CBOP
C    !ROUTINE: RECIP_HFACS_MACROS.h
C    !INTERFACE:
C    include RECIP_HFACS_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_HFACS_MACROS.h                                      
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP







C
CBOP
C    !ROUTINE: RECIP_HFACW_MACROS.h
C    !INTERFACE:
C    include RECIP_HFACW_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RECIP_HFACW_MACROS.h                                      
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP







C
CBOP
C    !ROUTINE: XC_MACROS.h
C    !INTERFACE:
C    include XC_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | XC_MACROS.h                                               
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: YC_MACROS.h
C    !INTERFACE:
C    include YC_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | YC_MACROS.h                                               
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: RA_MACROS.h
C    !INTERFACE:
C    include RA_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RA_MACROS.h                                               
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP




C
CBOP
C    !ROUTINE: RAW_MACROS.h
C    !INTERFACE:
C    include RAW_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RAW_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP




C
CBOP
C    !ROUTINE: RAS_MACROS.h
C    !INTERFACE:
C    include RAS_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | RAS_MACROS.h                                              
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: MASKW_MACROS.h
C    !INTERFACE:
C    include MASKW_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | MASKW_MACROS.h                                            
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP






C
CBOP
C    !ROUTINE: MASKS_MACROS.h
C    !INTERFACE:
C    include MASKS_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | MASKS_MACROS.h                                            
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP






C
CBOP
C    !ROUTINE: TANPHIATU_MACROS.h
C    !INTERFACE:
C    include TANPHIATU_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | TANPHIATU_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: TANPHIATV_MACROS.h
C    !INTERFACE:
C    include TANPHIATV_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | TANPHIATV_MACROS.h                                        
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C
CBOP
C    !ROUTINE: FCORI_MACROS.h
C    !INTERFACE:
C    include FCORI_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | FCORI_MACROS.h                                            
C     *==========================================================*
C     | These macros are used to reduce memory requirement and/or 
C     | memory references when variables are fixed along a given  
C     | axis or axes.                                             
C     *==========================================================*
C     \ev
CEOP





C--   COMMON /GRID_RL/ RL valued grid defining variables.
C     deepFacC  :: deep-model grid factor (fct of vertical only) for dx,dy
C     deepFacF     at level-center (deepFacC)  and level interface (deepFacF)
C     deepFac2C :: deep-model grid factor (fct of vertical only) for area dx*dy
C     deepFac2F    at level-center (deepFac2C) and level interface (deepFac2F)
C     gravitySign :: indicates the direction of gravity relative to R direction
C                   (= -1 for R=Z (Z increases upward, -gravity direction  )
C                   (= +1 for R=P (P increases downward, +gravity direction)
C     rkSign     :: Vertical coordinate to vertical index orientation.
C                   ( +1 same orientation, -1 opposite orientation )
C     globalArea :: Domain Integrated horizontal Area [m2]
      COMMON /GRID_RL/
     &  cosFacU, cosFacV, sqCosFacU, sqCosFacV,
     &  deepFacC, deepFac2C, recip_deepFacC, recip_deepFac2C,
     &  deepFacF, deepFac2F, recip_deepFacF, recip_deepFac2F,
     &  gravitySign, rkSign, globalArea
      Real*8 cosFacU        (1-OLy:sNy+OLy,nSx,nSy)
      Real*8 cosFacV        (1-OLy:sNy+OLy,nSx,nSy)
      Real*8 sqCosFacU      (1-OLy:sNy+OLy,nSx,nSy)
      Real*8 sqCosFacV      (1-OLy:sNy+OLy,nSx,nSy)
      Real*8 deepFacC       (Nr)
      Real*8 deepFac2C      (Nr)
      Real*8 deepFacF       (Nr+1)
      Real*8 deepFac2F      (Nr+1)
      Real*8 recip_deepFacC (Nr)
      Real*8 recip_deepFac2C(Nr)
      Real*8 recip_deepFacF (Nr+1)
      Real*8 recip_deepFac2F(Nr+1)
      Real*8 gravitySign
      Real*8 rkSign
      Real*8 globalArea

C--   COMMON /GRID_RS/ RS valued grid defining variables.
C     dxC     :: Cell center separation in X across western cell wall (m)
C     dxG     :: Cell face separation in X along southern cell wall (m)
C     dxF     :: Cell face separation in X thru cell center (m)
C     dxV     :: V-point separation in X across south-west corner of cell (m)
C     dyC     :: Cell center separation in Y across southern cell wall (m)
C     dyG     :: Cell face separation in Y along western cell wall (m)
C     dyF     :: Cell face separation in Y thru cell center (m)
C     dyU     :: U-point separation in Y across south-west corner of cell (m)
C     drC     :: Cell center separation along Z axis ( units of r ).
C     drF     :: Cell face separation along Z axis ( units of r ).
C     R_low   :: base of fluid in r_unit (Depth(m) / Pressure(Pa) at top Atmos.)
C     rLowW   :: base of fluid column in r_unit at Western  edge location.
C     rLowS   :: base of fluid column in r_unit at Southern edge location.
C     Ro_surf :: surface reference (at rest) position, r_unit.
C     rSurfW  :: surface reference position at Western  edge location [r_unit].
C     rSurfS  :: surface reference position at Southern edge location [r_unit].
C     hFac    :: Fraction of cell in vertical which is open i.e how
C              "lopped" a cell is (dimensionless scale factor).
C              Note: The code needs terms like MIN(hFac,hFac(I-1))
C                    On some platforms it may be better to precompute
C                    hFacW, hFacS, ... here than do MIN on the fly.
C     maskInC :: Cell Center 2-D Interior mask (i.e., zero beyond OB)
C     maskInW :: West  face 2-D Interior mask (i.e., zero on and beyond OB)
C     maskInS :: South face 2-D Interior mask (i.e., zero on and beyond OB)
C     maskC   :: cell Center land mask
C     maskW   :: West face land mask
C     maskS   :: South face land mask
C     recip_dxC   :: Reciprocal of dxC
C     recip_dxG   :: Reciprocal of dxG
C     recip_dxF   :: Reciprocal of dxF
C     recip_dxV   :: Reciprocal of dxV
C     recip_dyC   :: Reciprocal of dxC
C     recip_dyG   :: Reciprocal of dyG
C     recip_dyF   :: Reciprocal of dyF
C     recip_dyU   :: Reciprocal of dyU
C     recip_drC   :: Reciprocal of drC
C     recip_drF   :: Reciprocal of drF
C     recip_Rcol  :: Inverse of cell center column thickness (1/r_unit)
C     recip_hFacC :: Inverse of cell open-depth f[X,Y,Z] ( dimensionless ).
C     recip_hFacW    rhFacC center, rhFacW west, rhFacS south.
C     recip_hFacS    Note: This is precomputed here because it involves division.
C     xC        :: X-coordinate of cell center f[X,Y]. The units of xc, yc
C                  depend on the grid. They are not used in differencing or
C                  averaging but are just a convient quantity for I/O,
C                  diagnostics etc.. As such xc is in m for cartesian
C                  coordinates but degrees for spherical polar.
C     yC        :: Y-coordinate of center of cell f[X,Y].
C     yG        :: Y-coordinate of corner of cell ( c-grid vorticity point) f[X,Y].
C     rA        :: R-face are f[X,Y] ( m^2 ).
C                  Note: In a cartesian framework rA is simply dx*dy,
C                      however we use rA to allow for non-globally
C                      orthogonal coordinate frames (with appropriate
C                      metric terms).
C     rC        :: R-coordinate of center of cell f[Z] (units of r).
C     rF        :: R-coordinate of face of cell f[Z] (units of r).
C - *HybSigm* - :: Hybrid-Sigma vert. Coord coefficients
C     aHybSigmF    at level-interface (*HybSigmF) and level-center (*HybSigmC)
C     aHybSigmC    aHybSigm* = constant r part, bHybSigm* = sigma part, such as
C     bHybSigmF    r(ij,k,t) = rLow(ij) + aHybSigm(k)*[rF(1)-rF(Nr+1)]
C     bHybSigmC              + bHybSigm(k)*[eta(ij,t)+Ro_surf(ij) - rLow(ij)]
C     dAHybSigF :: vertical increment of Hybrid-Sigma coefficient: constant r part,
C     dAHybSigC    between interface (dAHybSigF) and between center (dAHybSigC)
C     dBHybSigF :: vertical increment of Hybrid-Sigma coefficient: sigma part,
C     dBHybSigC    between interface (dBHybSigF) and between center (dBHybSigC)
C     tanPhiAtU :: tan of the latitude at U point. Used for spherical polar
C                  metric term in U equation.
C     tanPhiAtV :: tan of the latitude at V point. Used for spherical polar
C                  metric term in V equation.
C     angleCosC :: cosine of grid orientation angle relative to Geographic direction
C               at cell center: alpha=(Eastward_dir,grid_uVel_dir)=(North_d,vVel_d)
C     angleSinC :: sine   of grid orientation angle relative to Geographic direction
C               at cell center: alpha=(Eastward_dir,grid_uVel_dir)=(North_d,vVel_d)
C     u2zonDir  :: cosine of grid orientation angle at U point location
C     v2zonDir  :: minus sine of  orientation angle at V point location
C     fCori     :: Coriolis parameter at grid Center point
C     fCoriG    :: Coriolis parameter at grid Corner point
C     fCoriCos  :: Coriolis Cos(phi) parameter at grid Center point (for NH)

      COMMON /GRID_RS/
     &  dxC,dxF,dxG,dxV,dyC,dyF,dyG,dyU,
     &  R_low, rLowW, rLowS,
     &  Ro_surf, rSurfW, rSurfS,
     &  hFacC, hFacW, hFacS,
     &  recip_dxC,recip_dxF,recip_dxG,recip_dxV,
     &  recip_dyC,recip_dyF,recip_dyG,recip_dyU,
     &  recip_Rcol,
     &  recip_hFacC,recip_hFacW,recip_hFacS,
     &  xC,yC,rA,rAw,rAs,rAz,xG,yG,
     &  maskInC, maskInW, maskInS,
     &  maskC, maskW, maskS,
     &  recip_rA,recip_rAw,recip_rAs,recip_rAz,
     &  drC, drF, recip_drC, recip_drF, rC, rF,
     &  aHybSigmF, bHybSigmF, aHybSigmC, bHybSigmC,
     &  dAHybSigF, dBHybSigF, dBHybSigC, dAHybSigC,
     &  tanPhiAtU, tanPhiAtV,
     &  angleCosC, angleSinC, u2zonDir, v2zonDir,
     &  fCori, fCoriG, fCoriCos
      Real*8 dxC            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 dxF            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 dxG            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 dxV            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 dyC            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 dyF            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 dyG            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 dyU            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 R_low          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 rLowW          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 rLowS          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 Ro_surf        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 rSurfW         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 rSurfS         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 hFacC          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      Real*8 hFacW          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      Real*8 hFacS          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      Real*8 recip_dxC      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_dxF      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_dxG      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_dxV      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_dyC      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_dyF      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_dyG      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_dyU      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_Rcol     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_hFacC    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      Real*8 recip_hFacW    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      Real*8 recip_hFacS    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      Real*8 xC             (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 xG             (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 yC             (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 yG             (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 rA             (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 rAw            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 rAs            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 rAz            (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_rA       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_rAw      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_rAs      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 recip_rAz      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 maskInC        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 maskInW        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 maskInS        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 maskC          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      Real*8 maskW          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      Real*8 maskS          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      Real*8 drC            (Nr+1)
      Real*8 drF            (Nr)
      Real*8 recip_drC      (Nr+1)
      Real*8 recip_drF      (Nr)
      Real*8 rC             (Nr)
      Real*8 rF             (Nr+1)
      Real*8 aHybSigmF      (Nr+1)
      Real*8 bHybSigmF      (Nr+1)
      Real*8 aHybSigmC      (Nr)
      Real*8 bHybSigmC      (Nr)
      Real*8 dAHybSigF      (Nr)
      Real*8 dBHybSigF      (Nr)
      Real*8 dBHybSigC      (Nr+1)
      Real*8 dAHybSigC      (Nr+1)
      Real*8 tanPhiAtU      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 tanPhiAtV      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 angleCosC      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 angleSinC      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 u2zonDir       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 v2zonDir       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 fCori          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 fCoriG         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 fCoriCos       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


C--   COMMON /GRID_I/ INTEGER valued grid defining variables.
C     kSurfC  :: vertical index of the surface tracer cell
C     kSurfW  :: vertical index of the surface U point
C     kSurfS  :: vertical index of the surface V point
C     kLowC   :: index of the r-lowest "wet cell" (2D)
C IMPORTANT: kLowC = 0 and kSurfC,W,S = Nr+1 (or =Nr+2 on a thin-wall)
C            where the fluid column is empty (continent)
      COMMON /GRID_I/
     &  kSurfC, kSurfW, kSurfS,
     &  kLowC
      INTEGER kSurfC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER kSurfW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER kSurfS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER kLowC (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C $Header: /u/gcmpack/MITgcm/verification/offline_exf_seaice/code/SEAICE_SIZE.h,v 1.4 2014/06/25 10:54:27 mlosch Exp $
C $Name:  $


CBOP
C    !ROUTINE: SEAICE_SIZE.h
C    !INTERFACE:
C #include SEAICE_SIZE.h

C    !DESCRIPTION:
C Contains seaice array-size definition (number of tracers,categories).

C SItrMaxNum :: number of passive tracers to allocate
C nITD       :: number of seaice categories to allocate
CEOP

C-    Maximum Number of categories
      INTEGER nITD
C--
      PARAMETER (nITD = 1)

C-    Maximum Number of tracers
      INTEGER SItrMaxNum
      PARAMETER(SItrMaxNum = 3 )



CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
C     *==========================================================*
C     | SEAICE_PARAMS.h
C     | o Basic parameter header for sea ice model.
C     *==========================================================*

C--   COMMON /SEAICE_PARM_L/ Logical parameters of sea ice model.
C - dynamics:
C     SEAICEuseDYNAMICS :: If false, do not use dynamics;
C                          default is to use dynamics.
C     SEAICEuseFREEDRIFT :: If True use free drift velocity instead of EVP
C                           or LSR
C     SEAICEuseStrImpCpl:: If true use strongly implicit coupling formulation
C                          for LSR solver (Hutchings et al 2004,
C                          Ocean Modelling, eq.44)
C     SEAICEuseEVP      :: If true use elastic viscous plastic solver
C     SEAICEuseEVPstar  :: If true use modified elastic viscous plastic
C                          solver (EVP*) by Lemieux et al (2012)
C     SEAICEuseEVPrev   :: If true use "revisited" elastic viscous plastic
C                          solver following Bouillon et al. (2013), very similar
C                          to EVP*, but uses fewer implicit terms and drops
C                          one 1/e^2 in equations for sigma2 and sigma12
C     SEAICEuseEVPpickup :: Set to false in order to start EVP solver with
C                          non-EVP pickup files.  Default is true.
C                          Applied only if SEAICEuseEVP=.TRUE.
C     SEAICEuseMultiTileSolver :: in LSR, use full domain tri-diagonal solver
C     SEAICEuseLSR      :: If true, use default Picard solver with Line-
C                          Successive(-over)-Relaxation, can also be true
C                          if LSR is used as a preconditioner for the
C                          non-linear JFNK solver
C     SEAICEusePicardAsPrecon :: If true, allow SEAICEuseLSR = .TRUE. as a
C                          preconditioner for non-linear JFNK problem (def. = F)
C     SEAICEuseKrylov   :: If true, use matrix-free Krylov solver with Picard
C                          solver instead of LSR (default: false)
C     SEAICEuseJFNK     :: If true, use Jacobi-free Newton-Krylov solver
C                          instead of LSR (default: false)
C     SEAICEuseIMEX     :: use IMplicit/EXplicit scheme with JFNK
C     SEAICEuseTEM      :: to use the truncated ellipse method (Geiger et al.)
C                          1998) set this parameter to true, default is false
C     SEAICEuseMCS      :: to use the Mohr-Coulomb yield curve with a shear
C                          only flow rule (Ip et al 1991), set this parameter to
C                          true, default is false
C     SEAICEuseMCE      :: to use the Mohr-Coulomb yield curve with elliptical
C                          plastic potential (similarly to Hibler and Schulson
C                          2000 without the elliptical cap) set this parameter
C                          to true, default is false
C     SEAICEuseTD       :: to use the teardrop yield curve (Zhang and Rothrock,
C                          2005) set this parameter to true, default is false
C     SEAICEusePL       :: to use the parabolic lens yield curve (Zhang and
C                          Rothrock, 2005) set this parameter to true,
C                          default is false
C     SEAICEuseMEB      :: use Maxwell-Elasto-Brittle rheology (within
C                          implicity solver for momentum equations)
C     SEAICEupdateDamage:: update damage parameter (true for MEB,
C                          false otherwise)
C     SEAICEuseTilt     :: If true then include surface tilt term in dynamics
C     SEAICEuseMetricTerms :: use metric terms for dynamics solver
C                          (default = .true. )
C     SEAICE_no_slip    :: apply no slip boundary conditions to seaice velocity
C     SEAICE_2ndOrderBC :: apply 2nd order no slip boundary conditions (works
C                          only with EVP, JFNK or KRYLOV solver, default=F)
C     SEAICE_maskRHS    :: mask the RHS of the solver where there is no ice
C     SEAICE_clipVelocities :: clip velocities to +/- 40cm/s
C     SEAICEaddSnowMass :: in computing seaiceMass, add snow contribution
C                          default is .TRUE.
C     useHB87stressCoupling :: use an intergral over ice and ocean surface
C                          layer to define surface stresses on ocean
C                          following Hibler and Bryan (1987, JPO)
C     SEAICEupdateOceanStress :: If TRUE, update ocean surface stress
C                                accounting for seaice cover (default= T)
C     SEAICEuseBDF2     :: use 2nd-order backward difference approach
C                          for momentum equations as described in
C                          Lemieux et al. 2014, JCP
C                          so far only implemented for JFNK-solver
C     useHibler79IceStrength :: if true original ice strength parameterization
C                          other use Rothrock (1975) parameterization based
C                          on energetics and an ice thickness distribution
C                          (default = .true.)
C     SEAICEscaleSurfStress :: if TRUE, scale ice-ocean and ice-atmosphere
C                          stress on ice by concenration (AREA) following
C                          Connolley et al. (2004), JPO. (default = .TRUE.)
C     SEAICEsimpleRidging :: use Hibler(1979) ridging (default=.true.)
C     SEAICEuseLinRemapITD :: use linear remapping (Lipscomb et al. 2001)
C                             .TRUE. by default
C - advection:
C     SEAICEuseFluxForm :: use flux form for advection and diffusion
C                          of seaice
C     SEAICEadvHeff     :: turn on advection of effective thickness
C                          (default = .true.)
C     SEAICEadvArea     :: turn on advection of fraction area
C                          (default = .true.)
C     SEAICEadvSnow     :: turn on advection of snow (does not work with
C                          non-default Leap-frog scheme for advection)
C     SEAICEadvSalt     :: turn on advection of salt (does not work with
C                          non-default Leap-frog scheme for advection)
C     SEAICEmultiDimAdvection:: internal flag, set to false if any sea ice
C                          variable uses a non-multi-dimensional advection
C                          scheme
C     SEAICEmomAdvection:: turn on advection of momentum (default = .false.)
C     SEAICEhighOrderVorticity :: momentum advection parameters analogous to
C     SEAICEupwindVorticity    :: highOrderVorticity, upwindVorticity,
C     SEAICEuseAbsVorticity    :: useAbsVorticity, useJamartMomAdv for vector
C     SEAICEuseJamartMomAdv    :: invariant momentum in the ocean
C - thermodynamics:
C     usePW79thermodynamics :: use "0-layer" thermodynamics as described in
C                           Parkinson and Washington (1979) and Hibler (1979)
C     SEAICE_useMultDimSnow :: use same fixed pdf for snow as for
C                              multi-thickness-category ice (default=.TRUE.)
C     SEAICEuseFlooding  :: turn on scheme to convert submerged snow into ice
C     SEAICEheatConsFix  :: If true then fix ocn<->seaice advective heat flux.
C     useMaykutSatVapPoly :: use Maykut Polynomial for saturation vapor pressure
C                         instead of extended temp-range exponential law; def=F.
C     SEAICE_mcPheeStepFunc :: use step function (not linear tapering) in
C                           ocean-ice turbulent flux
C     SEAICE_doOpenWaterGrowth :: use open water heat flux directly to grow ice
C                           (when false cool ocean, and grow later if needed)
C     SEAICE_doOpenWaterMelt   :: use open water heat flux directly to melt ice
C                           (when false warm ocean, and melt later if needed)
C     SEAICE_growMeltByConv :: grow/melt according to convergence of turbulence
C                              and conduction, rather than in two steps (default)
C     SEAICE_salinityTracer    :: use SItracer to exchange and trace ocean
C                           salt in ice
C     SEAICE_ageTracer         :: use SItracer to trace the age of ice
C     SEAICErestoreUnderIce :: restore surface T/S also underneath ice
C                          ( default is false )
C - other (I/O, ...):
C     SEAICEwriteState  :: If true, write sea ice state to file;
C                          default is false.
C     SEAICE_tave_mdsio :: write TimeAverage output using MDSIO
C     SEAICE_dump_mdsio :: write snap-shot output   using MDSIO
C     SEAICE_mon_stdio  :: write monitor to std-outp
C     SEAICE_tave_mnc   :: write TimeAverage output using MNC
C     SEAICE_dump_mnc   :: write snap-shot output   using MNC
C     SEAICE_mon_mnc    :: write monitor to netcdf file
      LOGICAL
     &     SEAICEuseDYNAMICS, SEAICEuseFREEDRIFT, SEAICEuseStrImpCpl,
     &     SEAICEuseEVP, SEAICEuseEVPstar, SEAICEuseEVPrev,
     &     SEAICEuseEVPpickup,
     &     SEAICEuseMultiTileSolver,
     &     SEAICEuseLSR, SEAICEuseKrylov,
     &     SEAICEuseJFNK, SEAICEuseIMEX, SEAICEuseBDF2,SEAICEuseMEB,
     &     SEAICEupdateDamage,
     &     SEAICEusePicardAsPrecon,
     &     useHibler79IceStrength, SEAICEsimpleRidging,
     &     SEAICEuseLinRemapITD, SEAICEuseTD, SEAICEusePL,
     &     SEAICEuseTEM, SEAICEuseTilt, SEAICEuseMetricTerms,
     &     SEAICEuseMCS, SEAICEuseMCE,
     &     SEAICE_no_slip, SEAICE_2ndOrderBC,
     &     SEAICE_maskRHS, SEAICEscaleSurfStress,
     &     SEAICE_clipVelocities, SEAICEaddSnowMass,
     &     useHB87stressCoupling, SEAICEupdateOceanStress,
     &     SEAICEuseFluxForm, SEAICEadvHeff, SEAICEadvArea,
     &     SEAICEmultiDimAdvection,
     &     SEAICEadvSnow, SEAICEadvSalt, SEAICEmomAdvection,
     &     SEAICEhighOrderVorticity, SEAICEupwindVorticity,
     &     SEAICEuseAbsVorticity, SEAICEuseJamartMomAdv,
     &     usePW79thermodynamics,
     &     SEAICE_useMultDimSnow, SEAICEuseFlooding, SEAICEheatConsFix,
     &     useMaykutSatVapPoly, SEAICE_mcPheeStepFunc,
     &     SEAICE_doOpenWaterGrowth, SEAICE_doOpenWaterMelt,
     &     SEAICE_salinityTracer, SEAICE_ageTracer,
     &     SEAICErestoreUnderIce, SEAICE_growMeltByConv,
     &     SEAICEwriteState,
     &     SEAICE_tave_mdsio, SEAICE_dump_mdsio, SEAICE_mon_stdio,
     &     SEAICE_tave_mnc,   SEAICE_dump_mnc,   SEAICE_mon_mnc
      COMMON /SEAICE_PARM_L/
     &     SEAICEuseDYNAMICS, SEAICEuseFREEDRIFT, SEAICEuseStrImpCpl,
     &     SEAICEuseEVP, SEAICEuseEVPstar, SEAICEuseEVPrev,
     &     SEAICEuseEVPpickup,
     &     SEAICEuseMultiTileSolver,
     &     SEAICEuseLSR, SEAICEuseKrylov,
     &     SEAICEuseJFNK, SEAICEuseIMEX, SEAICEuseBDF2, SEAICEuseMEB,
     &     SEAICEupdateDamage,
     &     SEAICEusePicardAsPrecon,
     &     useHibler79IceStrength, SEAICEsimpleRidging,
     &     SEAICEuseLinRemapITD, SEAICEuseTD, SEAICEusePL,
     &     SEAICEuseTEM, SEAICEuseTilt, SEAICEuseMetricTerms,
     &     SEAICEuseMCS, SEAICEuseMCE,
     &     SEAICE_no_slip, SEAICE_2ndOrderBC,
     &     SEAICE_maskRHS, SEAICEscaleSurfStress,
     &     SEAICE_clipVelocities, SEAICEaddSnowMass,
     &     useHB87stressCoupling, SEAICEupdateOceanStress,
     &     SEAICEuseFluxForm, SEAICEadvHeff, SEAICEadvArea,
     &     SEAICEadvSnow, SEAICEadvSalt, SEAICEmomAdvection,
     &     SEAICEmultiDimAdvection,
     &     SEAICEhighOrderVorticity, SEAICEupwindVorticity,
     &     SEAICEuseAbsVorticity, SEAICEuseJamartMomAdv,
     &     usePW79thermodynamics,
     &     SEAICE_useMultDimSnow, SEAICEuseFlooding, SEAICEheatConsFix,
     &     useMaykutSatVapPoly, SEAICE_mcPheeStepFunc,
     &     SEAICE_doOpenWaterGrowth, SEAICE_doOpenWaterMelt,
     &     SEAICE_salinityTracer, SEAICE_ageTracer,
     &     SEAICErestoreUnderIce, SEAICE_growMeltByConv,
     &     SEAICEwriteState,
     &     SEAICE_tave_mdsio, SEAICE_dump_mdsio, SEAICE_mon_stdio,
     &     SEAICE_tave_mnc,   SEAICE_dump_mnc,   SEAICE_mon_mnc

C--   COMMON /SEAICE_PARM_I/ Integer valued parameters of sea ice model.
C     IMAX_TICE         :: number of iterations for ice surface temp
C                          (default=10)
C     postSolvTempIter :: select flux calculation after surf. temp solver
C                         iteration
C                         0 = none, i.e., from last iter
C                         1 = use linearized approx (consistent with tsurf
C                             finding)
C                         2 = full non-lin form
C     SOLV_NCHECK         :: iteration interval for LSR-solver convergence test
C     SEAICEnonLinIterMax :: number of allowed non-linear solver iterations
C                            for implicit solvers (JFNK and Picard) (>= 2)
C     SEAICElinearIterMax :: number of allowed linear solver iterations for
C                            for implicit solvers (JFNK and Picard) C
C     SEAICEpreconNL_Iter :: number non-linear iterations in preconditioner
C     SEAICEpreconLinIter :: number linear iterations in preconditioner
C     SEAICEnEVPstarSteps :: number of evp*-steps
C     SEAICEmomStartBDF   :: number of previous u/vIce time levels available
C                          to start (or restart) BDF2 scheme.
C     SEAICE_JFNK_lsIter  :: number of Newton iterations after which the
C                            line search is started
C     SEAICE_JFNK_tolIter :: number of Newton iterations after which the
C                            the tolerance is relaxed again (default = 100)
C     SEAICE_OLx/y      :: overlaps for LSR-solver and for the
C                          LSR-preconditioner in JFNK and KRYLOV solver;
C                          for 0 < SEAICE_OLx/y 0 <= OLx/y-2 the LSR solver
C                          and preconditioner use a restricted additive
C                          Schwarz method (default = OLx/y-2).
C     LSR_mixIniGuess   :: control mixing of free-drift sol. into LSR initial
C                          guess
C                       :: =0 : nothing; =1 : no mix, but print free-drift
C                          resid.;
C                       :: =2,4 : mix with (1/local-residual)^2,4 factor
C     SEAICEpresPow0    :: HEFF exponent for ice strength below SEAICEpresH0
C     SEAICEpresPow1    :: HEFF exponent for ice strength above SEAICEpresH0
C     rigding parameters (only active when SEAICE_ITD is defined)
C     SEAICEpartFunc    :: =0 use Thorndyke et al (1975) participation function
C                       :: =1 use Lipscomb et al (2007) participation function
C     SEAICEredistFunc  :: =0 assume ridged ice is uniformly distributed
C                             (Hibler, 1980)
C                          =1 Following Lipscomb et al. (2007), ridged ice is
C                             distributed following an exponentially
C                             decaying function
C     SEAICEridgingIterMax :: maximum number of ridging iterations
C     end ridging parameters
C     SEAICEselectKEscheme   :: momentum advection parameters analogous
C     SEAICEselectVortScheme :: to selectKEscheme and selectVortScheme
C     SEAICEadvScheme   :: sets the advection scheme for thickness and area
C                          (default = 77)
C     SEAICEadvSchArea  :: sets the advection scheme for area
C     SEAICEadvSchHeff  :: sets the advection scheme for effective thickness
C                         (=volume), snow thickness, and salt if available
C     SEAICEadvSchSnow  :: sets the advection scheme for snow on sea-ice
C     SEAICEadvSchSalt  :: sets the advection scheme for sea ice salinity
C     SEAICEadvSchSnow  :: sets the advection scheme for snow on sea-ice
C     SEAICE_areaLossFormula :: selects formula for ice cover loss from melt
C                        :: 1=from all but only melt conributions by ATM and OCN
C                        :: 2=from net melt-growth>0 by ATM and OCN
C                        :: 3=from predicted melt by ATM
C     SEAICE_areaGainFormula :: selects formula for ice cover gain from open
C                               water growth
C                        :: 1=from growth by ATM
C                        :: 2=from predicted growth by ATM
C     SEAICEetaZmethod   :: determines how shear-viscosity eta is computed at
C                           Z-points
C                           0=simple averaging from C-points (default and old)
C                           3=weighted averaging of squares of strain rates
C                             (recommended for energy conservation)
C     SEAICE_multDim     :: number of ice categories
C     SEAICE_debugPointI :: I,J index for seaice-specific debuggin
C     SEAICE_debugPointJ
C
      INTEGER IMAX_TICE, postSolvTempIter
      INTEGER SOLV_NCHECK
      INTEGER SEAICEnonLinIterMax, SEAICElinearIterMax
      INTEGER SEAICEpreconLinIter, SEAICEpreconNL_Iter
      INTEGER LSR_mixIniGuess
      INTEGER SEAICEnEVPstarSteps
      INTEGER SEAICEmomStartBDF
      INTEGER SEAICE_JFNK_lsIter, SEAICE_JFNK_tolIter
      INTEGER SEAICE_OLx, SEAICE_OLy
      INTEGER SEAICEselectKEscheme, SEAICEselectVortScheme
      INTEGER SEAICEadvScheme
      INTEGER SEAICEadvSchArea
      INTEGER SEAICEadvSchHeff
      INTEGER SEAICEadvSchSnow
      INTEGER SEAICEadvSchSalt
      INTEGER SEAICEadjMODE
      INTEGER SEAICE_areaLossFormula
      INTEGER SEAICE_areaGainFormula
      INTEGER SEAICEetaZmethod
      INTEGER SEAICE_multDim
      INTEGER SEAICE_debugPointI
      INTEGER SEAICE_debugPointJ
      INTEGER SEAICEpresPow0, SEAICEpresPow1
      INTEGER SEAICEpartFunc, SEAICEredistFunc
      INTEGER SEAICEridgingIterMax
      COMMON /SEAICE_PARM_I/
     &     IMAX_TICE, postSolvTempIter, SOLV_NCHECK,
     &     SEAICEnonLinIterMax, SEAICElinearIterMax,
     &     SEAICEpreconLinIter, SEAICEpreconNL_Iter,
     &     LSR_mixIniGuess,
     &     SEAICEnEVPstarSteps,
     &     SEAICEmomStartBDF,
     &     SEAICE_JFNK_lsIter, SEAICE_OLx, SEAICE_OLy,
     &     SEAICE_JFNK_tolIter,
     &     SEAICEpresPow0, SEAICEpresPow1,
     &     SEAICEpartFunc, SEAICEredistFunc, SEAICEridgingIterMax,
     &     SEAICEselectKEscheme, SEAICEselectVortScheme,
     &     SEAICEadvScheme,
     &     SEAICEadvSchArea,
     &     SEAICEadvSchHeff,
     &     SEAICEadvSchSnow,
     &     SEAICEadvSchSalt,
     &     SEAICEadjMODE,
     &     SEAICE_areaLossFormula,
     &     SEAICE_areaGainFormula,
     &     SEAICE_multDim,
     &     SEAICEetaZmethod,
     &     SEAICE_debugPointI,
     &     SEAICE_debugPointJ

C--   COMMON /SEAICE_PARM_C/ Character valued sea ice model parameters.
C     AreaFile          :: File containing initial sea-ice concentration
C     HsnowFile         :: File containing initial snow thickness
C     HsaltFile         :: File containing initial sea ice salt content
C     HeffFile          :: File containing initial sea-ice thickness
C     uIceFile          :: File containing initial sea-ice U comp. velocity
C     vIceFile          :: File containing initial sea-ice V comp. velocity
C        !!! NOTE !!! Initial sea-ice thickness can also be set using
C        SEAICE_initialHEFF below.  But a constant initial condition
C        can mean large artificial fluxes of heat and freshwater in
C        the surface layer during the first model time step.
C
      CHARACTER*(MAX_LEN_FNAM) AreaFile
      CHARACTER*(MAX_LEN_FNAM) HsnowFile
      CHARACTER*(MAX_LEN_FNAM) HsaltFile
      CHARACTER*(MAX_LEN_FNAM) HeffFile
      CHARACTER*(MAX_LEN_FNAM) uIceFile
      CHARACTER*(MAX_LEN_FNAM) vIceFile
      COMMON /SEAICE_PARM_C/
     &   AreaFile, HsnowFile, HsaltFile, HeffFile,
     &   uIceFile, vIceFile

C--   COMMON /SEAICE_PARM_RL/ Real valued parameters of sea ice model.
C     SEAICE_deltaTtherm :: Seaice timestep for thermodynamic equations (s)
C     SEAICE_deltaTdyn   :: Seaice timestep for dynamic solver          (s)
C     SEAICE_LSRrelaxU/V :: relaxation parameter for LSR-solver: U/V-component
C     SEAICE_deltaTevp   :: Seaice timestep for EVP solver              (s)
C     SEAICE_elasticParm :: parameter that sets relaxation timescale
C                           tau = SEAICE_elasticParm * SEAICE_deltaTdyn
C     SEAICE_evpTauRelax :: relaxation timescale tau                    (s)
C     SEAICE_evpDampC    :: evp damping constant (Hunke,JCP,2001)       (kg/m^2)
C     SEAICE_evpAlpha    :: dimensionless parameter 2*evpTauRelax/deltaTevp
C     SEAICE_evpBeta     :: dimensionless parameter deltaTdyn/deltaTevp
C     SEAICEaEVPcoeff    :: main coefficent for adaptive EVP (largest
C                           stabilized frequency)
C     SEAICEaEVPcStar    :: multiple of stabilty factor: alpha*beta=cstar*gamma
C     SEAICEaEVPalphaMin :: lower limit of alpha and beta, regularisation
C                           to prevent singularities of system matrix,
C                           e.g. when ice concentration is too low.
C     SEAICEnonLinTol    :: non-linear tolerance parameter for implicit solvers
C     JFNKgamma_lin_min/max :: tolerance parameters for linear JFNK solver
C     JFNKres_t          :: tolerance parameter for FGMRES residual
C     JFNKres_tFac       :: if set, JFNKres_t=JFNKres_tFac*(initial residual)
C     SEAICE_JFNKepsilon :: step size for the FD-gradient in s/r seaice_jacvec
C     SEAICE_JFNKphi     :: [0,1] parameter for inexact Newton Method (def = 1)
C     SEAICE_JFNKalpha   :: (1,2] parameter for inexact Newton Method (def = 1)
C     SEAICE_zetaMaxFac  :: factor determining the maximum viscosity    (s)
C                          (default = 5.e+12/2.e4 = 2.5e8)
C     SEAICE_zetaMin     :: lower bound for viscosity (default = 0)    (N s/m^2)
C     SEAICEpresH0       :: HEFF threshold for ice strength            (m)
C     SEAICE_monFreq     :: SEAICE monitor frequency.                   (s)
C     SEAICE_dumpFreq    :: SEAICE dump frequency.                      (s)
C     SEAICE_taveFreq    :: SEAICE time-averaging frequency.            (s)
C     SEAICE_initialHEFF :: initial sea-ice thickness                   (m)
C     SEAICE_rhoAir      :: density of air                              (kg/m^3)
C     SEAICE_rhoIce      :: density of sea ice                          (kg/m^3)
C     SEAICE_rhoSnow     :: density of snow                             (kg/m^3)
C     ICE2WATR           :: ratio of sea ice density to water density
C     SEAICE_cpAir       :: specific heat of air                        (J/kg/K)
C
C     SEAICEpoissonRatio :: Maxwell-Elasto-Brittle parameter (0.3)
C     SEAICEdamageMin    :: lower cut off for regularization (1e-16)
C     SEAICEviscosity    :: SEAICE_strength * 1e7
C     SEAICEcohesion     :: constant cohesion parameter: 25 kN/m^2
C     SEAICEdamageParm   :: parameter alpha that determines ratio of damaged
C                           elasticity modululs undamaged appararent viscosity
C                           > 1 (default = 4.0)
C     SEAICEintFrictCoeff:: used to compute Mohr-Coulomb slope in principle
C                           stress space (default = 0.7)
C     SEAICEmohrCoulombSlope :: slope in stress invariant space,
C                           e.g. cos(45) = 1/sqrt(2)
C     SEAICEhealingTime  :: 1.D5
C     SEAICEdamageTime   :: 2.D0
C
C     OCEAN_drag         :: unitless air-ocean drag coefficient (default 0.001)
C     SEAICE_drag        :: unitless air-ice drag coefficient   (default 0.001)
C     SEAICE_waterDrag   :: unitless water-ice drag coefficient (default 0.0055)
C     SEAICEdWatMin      :: minimum linear water-ice drag applied to DWATN
C                           (default 0.25 m/s)
C
C     SEAICE_dryIceAlb   :: winter albedo
C     SEAICE_wetIceAlb   :: summer albedo
C     SEAICE_drySnowAlb  :: dry snow albedo
C     SEAICE_wetSnowAlb  :: wet snow albedo
C     HO                 :: AKA "lead closing parameter", demarcation thickness
C                           between thin and thick ice. Alternatively, HO (in
C                           meters) can be interpreted as the thickness of ice
C                           formed in open water.
C                           HO is a key ice-growth parameter that determines
C                           the partition between vertical and lateral growth.
C                           The default is 0.5m, increasing this value leads
C                           slower formation of a closed ice cover and thus to
C                           more ice (and thicker) ice, decreasing to faster
C                           formation of a closed ice cover (leads are closing
C                           faster) and thus less (thinner) ice.
C
C     SEAICE_drag_south       :: Southern Ocean SEAICE_drag
C     SEAICE_waterDrag_south  :: Southern Ocean SEAICE_waterDrag
C     SEAICE_dryIceAlb_south  :: Southern Ocean SEAICE_dryIceAlb
C     SEAICE_wetIceAlb_south  :: Southern Ocean SEAICE_wetIceAlb
C     SEAICE_drySnowAlb_south :: Southern Ocean SEAICE_drySnowAlb
C     SEAICE_wetSnowAlb_south :: Southern Ocean SEAICE_wetSnowAlb
C     HO_south                :: Southern Ocean HO
C
C     Parameters for basal drag of grounded ice following
C     Lemieux et al. (2015), doi:10.1002/2014JC010678
C     SEAICE_cBasalStar (default = SEAICE_cStar)
C     SEAICEbasalDragU0 (default = 5e-5)
C     SEAICEbasalDragK1 (default = 8)
C     SEAICEbasalDragK2  :: if > 0, turns on basal drag
C                           (default = 0, Lemieux suggests 15)
C
C     SEAICE_wetAlbTemp  :: Temp (deg.C) above which wet-albedo values are used
C     SEAICE_waterAlbedo :: water albedo
C     SEAICE_strength    :: sea-ice strength Pstar
C     SEAICE_cStar       :: sea-ice strength paramter C* (def: 20)
C     SEAICE_tensilFac   :: sea-ice tensile strength factor, values in [0,1]
C     SEAICE_tensilDepth :: crtical depth for sea-ice tensile strength (def 0.)
C     SEAICEpressReplFac :: interpolator between PRESS0 and regularized PRESS
C                           1. (default): pure pressure replace method (PRESS)
C                           0.          : pure Hibler (1979) method (PRESS0)
C     SEAICE_eccen       :: sea-ice eccentricity of the elliptical yield curve
C     SEAICE_eccfr       :: sea-ice eccentricity of the elliptical flow rule
C     SEAICE_lhFusion    :: latent heat of fusion for ice and snow (J/kg)
C     SEAICE_lhEvap      :: latent heat of evaporation for water (J/kg)
C     SEAICE_dalton      :: Dalton number (= sensible heat transfer coefficient)
C     SEAICE_iceConduct  :: sea-ice conductivity
C     SEAICE_snowConduct :: snow conductivity
C     SEAICE_emissivity  :: longwave ocean-surface emissivity (-)
C     SEAICE_ice_emiss   :: longwave ice-surface emissivity (-)
C     SEAICE_snow_emiss  :: longwave snow-surface emissivity (-)
C     SEAICE_boltzmann   :: Stefan-Boltzman constant (not a run time parameter)
C     SEAICE_snowThick   :: cutoff snow thickness (for snow-albedo)
C     SEAICE_shortwave   :: ice penetration shortwave radiation factor
C     SEAICE_saltFrac    :: salinity of newly formed seaice defined as a
C                           fraction of the ocean surface salinity at the time
C                           of freezing
C     SEAICE_salt0       :: prescribed salinity of seaice (in g/kg).
C     facOpenGrow        :: 0./1. version of logical SEAICE_doOpenWaterGrowth
C     facOpenMelt        :: 0./1. version of logical SEAICE_doOpenWaterMelt
C     SEAICE_mcPheePiston:: ocean-ice turbulent flux "piston velocity" (m/s)
C                           that sets melt efficiency.
C     SEAICE_mcPheeTaper :: tapering down of turbulent flux term with ice
C                           concentration. The 100% cover turb. flux is
C                           multiplied by 1.-SEAICE_mcPheeTaper
C     SEAICE_frazilFrac  :: Fraction of surface level negative heat content
C                           anomalies (relative to the local freezing point)
C                           may contribute as frazil over one time step.
C     SEAICE_tempFrz0    :: sea water freezing point is
C     SEAICE_dTempFrz_dS :: tempFrz = SEAICE_tempFrz0 + salt*SEAICE_dTempFrz_dS
C     SEAICE_PDF         :: prescribed sea-ice distribution within grid box
C     SEAICEstressFactor :: factor by which ice affects wind stress (default=1)
C     LSR_ERROR          :: sets accuracy of LSR solver
C     DIFF1              :: parameter used in advect.F
C     SEAICEtdMU         :: slope parameter for the teardrop and parabolic lens
C                           yield curves
C     SEAICE_deltaMin    :: small number used to reduce singularities of Delta
C     SEAICE_area_max    :: usually set to 1. Seeting areaMax below 1 specifies
C                           the minimun amount of leads (1-areaMax) in the
C                           ice pack.
C     SEAICE_area_floor  :: usually set to 1x10^-5. Specifies a minimun
C                           ice fraction in the ice pack.
C     SEAICE_area_reg    :: usually set to 1x10^-5. Specifies a minimun
C                           ice fraction for the purposes of regularization
C     SEAICE_hice_reg    :: usually set to 5 cm. Specifies a minimun
C                           ice thickness for the purposes of regularization
C     SEAICEdiffKhArea   :: sets the diffusivity for area (m^2/s)
C     SEAICEdiffKhHeff   :: sets the diffusivity for effective thickness (m^2/s)
C     SEAICEdiffKhSnow   :: sets the diffusivity for snow on sea-ice (m^2/s)
C     SEAICEdiffKhSalt   :: sets the diffusivity for sea ice salinity (m^2/s)
C     SEAICE_airTurnAngle   :: turning angles of air-ice interfacial stress
C     SEAICE_waterTurnAngle :: and ice-water interfacial stress (in degrees)
C     SEAICE_tauAreaObsRelax :: Timescale of relaxation to observed
C                               sea ice concentration (s), default=unset
C     ridging parameters (Lipscomb et al, 2007, Bitz et al. 2001):
C     SEAICE_cf       :: ratio of total energy sinks to gravitational sink
C                        (scales ice strength, suggested values: 2 to 17)
C     SEAICEgStar     :: maximum ice concentration that participates in ridging
C     SEAICEhStar     :: empirical thickness (ridging parameter)
C     SEAICEaStar     :: ice concentration parameter similar to gStar for
C                        exponential distribution (Lipscomb et al 2007)
C     SEAICEshearParm :: <=1 reduces amount of energy lost to ridge building
C     SEAICEmuRidging :: tuning parameter similar to hStar for Lipcomb et al
C                        (2007)-scheme
C     SEAICEmaxRaft   :: regularization parameter (default=1)
C     SEAICEsnowFracRidge :: fraction of snow that remains on ridged
C     SINegFac        :: SIADV over/undershoot factor in FW/Adjoint
C     SEAICEmcMu      :: parameter for MC yield curve for useMCE, useMCS and
C                        useTEM options, default is one
C
      Real*8 SEAICE_deltaTtherm, SEAICE_deltaTdyn, SEAICE_deltaTevp
      Real*8 SEAICE_LSRrelaxU, SEAICE_LSRrelaxV
      Real*8 SEAICE_monFreq, SEAICE_dumpFreq, SEAICE_taveFreq
      Real*8 SEAICE_initialHEFF
      Real*8 SEAICE_rhoAir, SEAICE_rhoIce, SEAICE_rhoSnow, ICE2WATR
      Real*8 SEAICE_cpAir
      Real*8 SEAICE_drag, SEAICE_waterDrag, SEAICEdWatMin
      Real*8 SEAICE_dryIceAlb, SEAICE_wetIceAlb
      Real*8 SEAICE_drySnowAlb, SEAICE_wetSnowAlb, HO
      Real*8 SEAICE_drag_south, SEAICE_waterDrag_south
      Real*8 SEAICE_dryIceAlb_south, SEAICE_wetIceAlb_south
      Real*8 SEAICE_drySnowAlb_south, SEAICE_wetSnowAlb_south, HO_south
      Real*8 SEAICE_cBasalStar, SEAICEbasalDragU0
      Real*8 SEAICEbasalDragK1, SEAICEbasalDragK2
      Real*8 SEAICE_wetAlbTemp, SEAICE_waterAlbedo
      Real*8 SEAICE_strength, SEAICE_cStar, SEAICEpressReplFac
      Real*8 SEAICE_tensilFac, SEAICE_tensilDepth
      Real*8 SEAICE_eccen, SEAICE_eccfr
      Real*8 SEAICEmcMu, SEAICEtdMU
      Real*8 SEAICE_lhFusion, SEAICE_lhEvap
      Real*8 SEAICE_dalton
      Real*8 SEAICE_iceConduct, SEAICE_snowConduct
      Real*8 SEAICE_emissivity, SEAICE_ice_emiss, SEAICE_snow_emiss
      Real*8 SEAICE_boltzmann
      Real*8 SEAICE_snowThick, SEAICE_shortwave
      Real*8 SEAICE_saltFrac, SEAICE_salt0, SEAICEstressFactor
      Real*8 SEAICE_mcPheeTaper, SEAICE_mcPheePiston
      Real*8 SEAICE_frazilFrac, SEAICE_availHeatFrac
      Real*8 facOpenGrow, facOpenMelt
      Real*8 SEAICE_tempFrz0, SEAICE_dTempFrz_dS
      Real*8 SEAICE_PDF(nITD)
      Real*8 OCEAN_drag, LSR_ERROR, DIFF1
      Real*8 SEAICEnonLinTol, JFNKres_t, JFNKres_tFac
      Real*8 JFNKgamma_lin_min, JFNKgamma_lin_max, SEAICE_JFNKepsilon
      Real*8 SEAICE_JFNKphi, SEAICE_JFNKalpha
      Real*8 SEAICE_deltaMin
      Real*8 SEAICE_area_reg, SEAICE_hice_reg
      Real*8 SEAICE_area_floor, SEAICE_area_max
      Real*8 SEAICE_airTurnAngle, SEAICE_waterTurnAngle
      Real*8 SEAICE_elasticParm, SEAICE_evpTauRelax
      Real*8 SEAICE_evpAlpha, SEAICE_evpBeta
      Real*8 SEAICE_evpDampC, SEAICE_zetaMin, SEAICE_zetaMaxFac
      Real*8 SEAICEaEVPcoeff, SEAICEaEVPcStar, SEAICEaEVPalphaMin
      Real*8 SEAICEpresH0
      Real*8 SEAICEdiffKhArea, SEAICEdiffKhHeff, SEAICEdiffKhSnow
      Real*8 SEAICEdiffKhSalt
      Real*8 SEAICE_tauAreaObsRelax
      Real*8 SEAICEgStar, SEAICEhStar, SEAICEaStar, SEAICEshearParm
      Real*8 SEAICEmuRidging, SEAICEmaxRaft, SEAICE_cf
      Real*8 SEAICEsnowFracRidge
      Real*8 SINegFac
      Real*8 SEAICEpoissonRatio, SEAICEviscosity, SEAICEcohesion
      Real*8 SEAICEdamageMin, SEAICEdamageParm
      Real*8 SEAICEintFrictCoeff, SEAICEmohrCoulombSlope
      Real*8 SEAICEhealingTime, SEAICEdamageTime

      COMMON /SEAICE_PARM_RL/
     &    SEAICE_deltaTtherm, SEAICE_deltaTdyn,
     &    SEAICE_LSRrelaxU, SEAICE_LSRrelaxV,
     &    SEAICE_deltaTevp, SEAICE_elasticParm, SEAICE_evpTauRelax,
     &    SEAICE_evpAlpha, SEAICE_evpBeta,
     &    SEAICEaEVPcoeff, SEAICEaEVPcStar, SEAICEaEVPalphaMin,
     &    SEAICE_evpDampC, SEAICE_zetaMin, SEAICE_zetaMaxFac,
     &    SEAICEpresH0,
     &    SEAICE_monFreq, SEAICE_dumpFreq, SEAICE_taveFreq,
     &    SEAICE_initialHEFF,
     &    SEAICE_rhoAir, SEAICE_rhoIce, SEAICE_rhoSnow, ICE2WATR,
     &    SEAICE_drag, SEAICE_waterDrag, SEAICEdWatMin,
     &    SEAICE_dryIceAlb, SEAICE_wetIceAlb,
     &    SEAICE_drySnowAlb, SEAICE_wetSnowAlb, HO,
     &    SEAICE_drag_south, SEAICE_waterDrag_south,
     &    SEAICE_dryIceAlb_south, SEAICE_wetIceAlb_south,
     &    SEAICE_drySnowAlb_south, SEAICE_wetSnowAlb_south, HO_south,
     &    SEAICE_cBasalStar, SEAICEbasalDragU0,
     &    SEAICEbasalDragK1, SEAICEbasalDragK2,
     &    SEAICE_wetAlbTemp, SEAICE_waterAlbedo,
     &    SEAICE_strength, SEAICE_cStar, SEAICE_eccen, SEAICE_eccfr,
     &    SEAICEtdMU, SEAICEmcMu,
     &    SEAICEpressReplFac, SEAICE_tensilFac, SEAICE_tensilDepth,
     &    SEAICE_lhFusion, SEAICE_lhEvap,
     &    SEAICE_dalton, SEAICE_cpAir,
     &    SEAICE_iceConduct, SEAICE_snowConduct,
     &    SEAICE_emissivity, SEAICE_ice_emiss, SEAICE_snow_emiss,
     &    SEAICE_boltzmann,
     &    SEAICE_snowThick, SEAICE_shortwave,
     &    SEAICE_saltFrac, SEAICE_salt0, SEAICEstressFactor,
     &    SEAICE_mcPheeTaper, SEAICE_mcPheePiston,
     &    SEAICE_frazilFrac, SEAICE_availHeatFrac,
     &    facOpenGrow, facOpenMelt,
     &    SEAICE_tempFrz0, SEAICE_dTempFrz_dS, SEAICE_PDF,
     &    OCEAN_drag, LSR_ERROR, DIFF1,
     &    SEAICEnonLinTol, JFNKres_t, JFNKres_tFac,
     &    JFNKgamma_lin_min, JFNKgamma_lin_max, SEAICE_JFNKepsilon,
     &    SEAICE_JFNKphi, SEAICE_JFNKalpha,
     &    SEAICE_deltaMin, SEAICE_area_reg, SEAICE_hice_reg,
     &    SEAICE_area_floor, SEAICE_area_max,
     &    SEAICEdiffKhArea, SEAICEdiffKhHeff, SEAICEdiffKhSnow,
     &    SEAICEdiffKhSalt, SEAICE_tauAreaObsRelax,
     &    SEAICE_airTurnAngle, SEAICE_waterTurnAngle,
     &    SEAICEgStar, SEAICEhStar, SEAICEaStar, SEAICEshearParm,
     &    SEAICEmuRidging, SEAICEmaxRaft, SEAICE_cf,
     &    SINegFac,
     &    SEAICEsnowFracRidge,
     &    SEAICEpoissonRatio, SEAICEviscosity, SEAICEcohesion,
     &    SEAICEdamageMin, SEAICEdamageParm,
     &    SEAICEintFrictCoeff, SEAICEmohrCoulombSlope,
     &    SEAICEhealingTime, SEAICEdamageTime

C--   COMMON /SEAICE_BOUND_RL/ Various bounding values
C     MIN_ATEMP         :: minimum air temperature   (deg C)
C     MIN_LWDOWN        :: minimum downward longwave (W/m^2)
C     MIN_TICE          :: minimum ice temperature   (deg C)
C     SEAICE_EPS        :: small number
C     SEAICE_EPS_SQ     :: small number square
C
      Real*8 MIN_ATEMP, MIN_LWDOWN, MIN_TICE
      Real*8 SEAICE_EPS, SEAICE_EPS_SQ
      COMMON /SEAICE_BOUND_RL/
     &     MIN_ATEMP, MIN_LWDOWN, MIN_TICE,
     &     SEAICE_EPS, SEAICE_EPS_SQ


C--   Constants used by sea-ice model
      Real*8         ZERO           , ONE           , TWO
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      Real*8         QUART            , HALF
      PARAMETER ( QUART = 0.25D0, HALF = 0.5D0 )
      Real*8 siEps
      PARAMETER ( siEps = 1.D-5 )

C--   Constants needed by McPhee formulas for turbulent ocean fluxes :
C        Stanton number (dimensionless), typical friction velocity
C        beneath sea ice (m/s), and tapering factor (dimensionless)
      Real*8 STANTON_NUMBER, USTAR_BASE, MCPHEE_TAPER_FAC
      PARAMETER ( MCPHEE_TAPER_FAC = 12.5D0 , STANTON_NUMBER =
     &            0.0056D0, USTAR_BASE = 0.0125D0 )

C--   identifiers for advected properties
      INTEGER GAD_HEFF,GAD_AREA,GAD_QICE1,GAD_QICE2,GAD_SNOW
      INTEGER GAD_SALT,GAD_SITR
      PARAMETER ( GAD_HEFF  = 1,
     &            GAD_AREA  = 2,
     &            GAD_SNOW  = 3,
     &            GAD_SALT  = 4,
     &            GAD_QICE1 = 5,
     &            GAD_QICE2 = 6,
     &            GAD_SITR  = 7)

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine Arguments ==
C     fld        :: In: input ice-field ; Out: laplacian of input field
C     DIFFA      :: grid length scale
C     iceMsk     :: Ocean/Land mask
C     myThid     :: my Thread Id. number
      Real*8 fld    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 DIFFA  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 iceMsk (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j,bi,bj :: Loop counters
      INTEGER i, j, bi, bj
      Real*8 DELTXX, DELTYY
      Real*8 tmpFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 dfx     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      Real*8 dfy     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        IF ( SEAICEuseFluxForm ) THEN
C--   Use flux form for MITgcm compliance, unfortunately changes results

          DO j=1-Oly,sNy+Oly
           DO i=1-Olx,sNx+Olx
            dfx(i,j) = 0.D0
            dfy(i,j) = 0.D0
           ENDDO
          ENDDO
C--   first compute fluxes across cell faces
          DO j=1,sNy+1
           DO i=1,sNx+1
            dfx(i,j) = dyG(i,j,bi,bj) * recip_dxC(i,j,bi,bj)
     &         * (fld(i,j,bi,bj)-fld(i-1,j,bi,bj))
     &         * cosFacU(j,bi,bj)
     &         * iceMsk(i,j,bi,bj)*iceMsk(i-1,j,bi,bj)
     &         * ( DIFFA(i,j,bi,bj)+DIFFA(i-1,j,bi,bj) )*HALF
            dfy(i,j) = dxG(i,j,bi,bj) * recip_dyC(i,j,bi,bj)
     &         * (fld(i,j,bi,bj)-fld(i,j-1,bi,bj))
     &         * iceMsk(i,j,bi,bj)*iceMsk(i,j-1,bi,bj)
     &         * ( DIFFA(i,j,bi,bj)+DIFFA(i,j-1,bi,bj) )*HALF
           ENDDO
          ENDDO
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            fld(i,j,bi,bj) = 0.D0
           ENDDO
          ENDDO
C--   compute Laplacian as flux divergence
          DO j=1,sNy
           DO i=1,sNx
            fld(i,j,bi,bj) = (
     &                         ( dfx(i+1,j) - dfx(i,j) )
     &                       + ( dfy(i,j+1) - dfy(i,j) )
     &                       ) * recip_rA(i,j,bi,bj)
           ENDDO
          ENDDO

        ELSE
C NOW DO DIFFUSION WITH NUIT CONVERSION

          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            tmpFld(i,j) = 0.0D0
           ENDDO
          ENDDO

          DO j=1,sNy
           DO i=1,sNx
            DELTXX = DIFFA(i,j,bi,bj)
     &             * recip_dxF(i,j,bi,bj)*recip_dxF(i,j,bi,bj)
            DELTYY = DIFFA(i,j,bi,bj)
     &             * recip_dyF(i,j,bi,bj)*recip_dyF(i,j,bi,bj)
     &             * recip_dxF(i,j,bi,bj)
            tmpFld(i,j) =
     &        DELTXX*(
     &                 (fld(i+1,j,bi,bj)-fld(i,  j,bi,bj))
     &                 *iceMsk(i+1,j,bi,bj)
     &                -(fld(i,  j,bi,bj)-fld(i-1,j,bi,bj))
     &                 *iceMsk(i-1,j,bi,bj)
     &               )
     &       +DELTYY*(
     &                 (fld(i,j+1,bi,bj)-fld(i,j,  bi,bj))
     &                 * dxG(i,j+1,bi,bj)*iceMsk(i,j+1,bi,bj)
     &                -(fld(i,j,  bi,bj)-fld(i,j-1,bi,bj))
     &                 * dxG(i,j,  bi,bj)*iceMsk(i,j-1,bi,bj)
     &               )
           ENDDO
          ENDDO
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            fld(i,j,bi,bj) = tmpFld(i,j)
           ENDDO
          ENDDO

C--  end flux-form / non flux-form
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END
