












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
C !ROUTINE: SEAICE_ITD_SUM

C !INTERFACE: ==========================================================
      SUBROUTINE SEAICE_ITD_SUM(
     I                  bi, bj, myTime, myIter, myThid )

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE SEAICE_ITD_SUM
C     | o sum ice area and volume over all ITD categories
C     |   and write into AREA and HEFF
C     |
C     | Torge Martin, Mar. 2012, torge@mit.edu
C     *===========================================================*
C \ev

C !USES: ===============================================================
      IMPLICIT NONE

C     === Global variables ===
C     nITD      :: number of seaice thickness categories
C     AREA      :: sea ice area   total
C     AREAITD   :: sea ice area   by category
C     HEFF      :: sea ice volume total
C     HEFFITD   :: sea ice volume by category
C     HSNOW     :: snow    volume total
C     HSNOWITD  :: snow    volume by category
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
c#include "PARAMS.h"
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
c#include "SEAICE_PARAMS.h"
CBOP
C !ROUTINE: SEAICE.h

C !DESCRIPTION: \bv
C     *==========================================================*
C     | SEAICE.h
C     | o Basic header for sea ice model.
C     |   Contains most sea ice field declarations.
C     *==========================================================*
C
C     UICE  :: zonal ice velocity in m/s at South-West B-grid
C              (or C-grid #ifdef ) U point
C              >0 from West to East
C     VICE  :: meridional ice velocity in m/s at South-West B-grid
C              (or C-grid #ifdef ) V point
C              >0 from South to North
C              note: the South-West B-grid U and V points are on
C                the lower, left-hand corner of each grid cell
C     AREA  :: fractional ice-covered area in m^2/m^2
C              at center of grid, i.e., tracer point
C              0 is no cover, 1 is 100% cover
C     HEFF  :: effective ice thickness in m
C              at center of grid, i.e., tracer point
C              note: for non-zero AREA, actual ice thickness is HEFF / AREA
C     HSNOW :: effective snow thickness in m
C              at center of grid, i.e., tracer point
C              note: for non-zero AREA, actual snow thickness is HSNOW / AREA
C     HSALT :: effective sea ice salinity in g/m^2
C              at center of grid, i.e., tracer point
C \ev
CEOP

C--   Grid variables for seaice
      COMMON/ARRAY/HEFFM, SIMaskU, SIMaskV
      Real*8 HEFFM      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     static masks (depend only on geometry)
      Real*8 SIMaskU    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 SIMaskV    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON/ARRAYC/ seaiceMaskU, seaiceMaskV
C     dynamic masks (depend on area)
      Real*8 seaiceMaskU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 seaiceMaskV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     k1/2AtZ :: coefficients at C and Z points
C     k1/2AtC    for metric terms in U/V ice equations.
      COMMON/ARRAYCMETRIC/  k1AtC, k1AtZ, k2AtC, k2AtZ
      Real*8 k1AtC      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 k1AtZ      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 k2AtC      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 k2AtZ      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C--   Dynamical variables
      COMMON/SEAICE_DYNVARS_1/AREA,HEFF,HSNOW,UICE,VICE
      Real*8 AREA       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 HEFF       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 HSNOW      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 UICE       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 VICE       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C     uIceC :: average of UICE between last two time steps
C     vIceC :: average of VICE between last two time steps
      COMMON/SEAICE_DYNVARS_3/
     &     ETA,etaZ,ZETA,zetaZ,PRESS, e11, e22, e12, deltaC,
     &     FORCEX,FORCEY,
     &     uIceNm1, vIceNm1
      Real*8 ETA        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 etaZ       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 ZETA       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 zetaZ      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     ice strength/pressure term
      Real*8 PRESS      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     strain rate tensor
      Real*8 e11        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 e22        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 e12        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     deformation rate tensor invariant, for viscous plastic sea ice =
C     sqrt[(e11**2+e22**2)*(1+1/e**2)+ 4./e**2*e12C**2 + 2*e11*e22*(1-1/e**2))
      Real*8 deltaC     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C
      Real*8 FORCEX     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 FORCEY     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 uIceNm1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 vIceNm1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


      COMMON/SEAICE_DYNVARS_CGRID/
     &     seaiceMassC, seaiceMassU, seaiceMassV
      Real*8 seaiceMassC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 seaiceMassU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 seaiceMassV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON/SEAICE_DYNVARS_4/
     &     DWATN, PRESS0, FORCEX0, FORCEY0, ZMAX, ZMIN, tensileStrFac
      Real*8 DWATN      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 PRESS0     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 FORCEX0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 FORCEY0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 ZMAX       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 ZMIN       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     factor k to compute the maximal tensile stress from k*PRESS0,
C     in analogy to the maximal compressive stress PRESS0
      Real*8 tensileStrFac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON/SEAICE_DAMAGE/ damage, damageNm1, cohesion, dCrit
      Real*8 damage     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 damageNm1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 dCrit      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 cohesion   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON/SEAICE_REG_NEG/d_HEFFbyNEG,d_HSNWbyNEG
C     The change of mean ice thickness due to out-of-bounds values following
C     sea ice dynamics and advection
      Real*8 d_HEFFbyNEG (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 d_HSNWbyNEG (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C

C     saltWtrIce contains m of salty ice melted (<0) or created (>0)
C     frWtrIce contains m of freshwater ice melted (<0) or created (>0)
C              that is, ice due to precipitation or snow
C     frWtrAtm contains freshwater flux from the atmosphere
      COMMON/ICEFLUX/ saltWtrIce, frWtrIce
      Real*8 saltWtrIce (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 frWtrIce   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C     TICES :: Seaice/snow surface temperature for each category
      COMMON/MULTICATEGORY/TICES
      Real*8 TICES      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nITD,nSx,nSy)


C
C     additional fields needed by the EVP solver and/or MEB rheology
C
C     seaice_sigma1  - sigma11+sigma22, defined at C-points
C     seaice_sigma2  - sigma11-sigma22, defined at C-points
C     seaice_sigma12 - off-diagonal term, defined at Z-points
      COMMON /SEAICE_EVP_FIELDS/
     &     seaice_sigma1, seaice_sigma2, seaice_sigma12
      Real*8 seaice_sigma1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 seaice_sigma2    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 seaice_sigma12   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     stressDivergenceX/Y - divergence of stress tensor
      COMMON /SEAICE_STRESSDIV/
     &     stressDivergenceX, stressDivergenceY
      Real*8 stressDivergenceX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8 stressDivergenceY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)



C     SWFracB :: fraction of surface Short-Wave radiation reaching
C                the bottom of ocean surface level
      Real*8 SWFracB
      COMMON /SEAICE_SW_R/
     &       SWFracB

C     diagnostics for the JFNK and Krylov solver
      INTEGER totalNewtonIters
      INTEGER totalNewtonFails
      INTEGER totalKrylovIters
      INTEGER totalKrylovFails
      INTEGER totalJFNKtimeSteps
      COMMON /SEAICE_SOLVER_I/
     &     totalNewtonIters, totalNewtonFails,
     &     totalKrylovIters, totalKrylovFails,
     &     totalJFNKtimeSteps
      INTEGER nVec
      PARAMETER ( nVec=2*sNx*sNy )
      Real*8 scalarProductMetric( nVec, 1, nSx, nSy )
      COMMON /SEAICE_KRYLOV_RL/ scalarProductMetric

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***

C !INPUT PARAMETERS: ===================================================
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      Real*8 myTime
      INTEGER bi, bj
      INTEGER myIter
      INTEGER myThid
CEOP

      RETURN
      END
