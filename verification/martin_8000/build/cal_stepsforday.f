












CBOP
C !ROUTINE: CAL_OPTIONS.h
C !INTERFACE:
C #include "CAL_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for Calendar (cal) package:
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


      subroutine cal_StepsForDay(
     I                            iday,
     O                            firststep,
     O                            laststep,
     O                            nsteps,
     I                            mythid
     &                          )

c     ==================================================================
c     SUBROUTINE cal_StepsForDay
c     ==================================================================
c
c     o Given the current day of the integration this routine returns
c       first, the last and the number of model steps the will have to
c       be performed.
c
c       This routine also checks consistency of variables quite
c       extensively.
c
c     started: Christian Eckert eckert@mit.edu  06-Apr-2000
c
c     changed:
c
c     ==================================================================
c     SUBROUTINE cal_StepsForDay
c     ==================================================================

      implicit none

c     == global variables ==

C     ==================================================================
C     HEADER calendar
C     ==================================================================
C
C     o This header file contains variables that are used by the
C       calendar tool. The calendar tool can be used in the ECCO
C       SEALION release of the MITgcmUV.
C
C     started: Christian Eckert eckert@mit.edu  30-Jun-1999
C     changed: Christian Eckert eckert@mit.edu  17-Dec-1999
C              - restructured the original version in order to have a
C                better interface to the MITgcmUV.
C
C     ==================================================================
C     HEADER calendar
C     ==================================================================

C   - Parameters of the numerical model:
C
C     modelStart       :: start time of the numerical model.
C     modelStartDate   :: start date of the numerical model.
C     modelEnd         :: end   time of the numerical model.
C     modelEndDate     :: end   date of the numerical model.
C     modelStep        :: timestep of the numerical model.
C     modelIntSteps    :: number of timestep that are to be performed.
C     modelIter0       :: the numerical models initial timestep number.
C     modelIterEnd     :: the models last timestep number.
C     modelStepsperday :: number of model time steps per day (<- removed).

C   - Parameters used by the calendar:
C
C     refDate          :: first day of the Gregorian Calendar.
C     nMonthYear       :: number months in a year.
C     nDayMonth        :: days per month depending on the year being a leap
C                         year or not. If the Model calendar is used a 360
C                         days year with 30 days months is used instead.
C     nDaysNoLeap      :: number of days in a usual year.
C     nDaysLeap        :: number of days in a leap year.
C     nMaxDayMonth     :: maximum number of days in a years month.
C     hoursPerDay      :: number of hours   in a calendars day.
C     minutesPerDay    :: number of minutes in a calendars day.
C     minutesPerHour   :: number of minutes in a calendars hour.
C     secondsPerDay    :: number of seconds in a calendars day.
C     secondsPerHour   :: number of seconds in a calendars hour.
C     secondsPerMinute :: number of seconds in a calendars minute.
C     cal_setStatus    :: status of calendar parms setting (0=none, 3=fully set)

      INTEGER nMonthYear
      PARAMETER ( nMonthYear = 12 )

      COMMON /CALENDAR_RL/
     &                modelStart,
     &                modelEnd,
     &                modelStep
      Real*8 modelStart
      Real*8 modelEnd
      Real*8 modelStep

      COMMON /CALENDAR_I/
     &               refDate,
     &               nDayMonth,
     &               nDaysNoLeap,
     &               nDaysLeap,
     &               nMaxDayMonth,
     &               hoursPerDay,
     &               minutesPerDay,
     &               minutesPerHour,
     &               secondsPerDay,
     &               secondsPerHour,
     &               secondsPerMinute,
     &               modelStartDate,
     &               modelEndDate,
     &               modelIter0,
     &               modelIterEnd,
     &               modelIntSteps,
     &               cal_setStatus,
     &               startdate_1,
     &               startdate_2

      INTEGER refDate(4)
      INTEGER nDayMonth(nMonthYear,2)
      INTEGER nDaysNoLeap
      INTEGER nDaysLeap
      INTEGER nMaxDayMonth
      INTEGER hoursPerDay
      INTEGER minutesPerDay
      INTEGER minutesPerHour
      INTEGER secondsPerDay
      INTEGER secondsPerHour
      INTEGER secondsPerMinute

      INTEGER modelStartDate(4)
      INTEGER modelEndDate(4)
      INTEGER modelIter0
      INTEGER modelIterEnd
      INTEGER modelIntSteps

      INTEGER cal_setStatus
      INTEGER startdate_1
      INTEGER startdate_2

C   calendarDumps :: When set, approximate months (30-31 days) and years (360-372 days)
C                    for parameters chkPtFreq, pChkPtFreq, taveFreq, SEAICE_taveFreq,
C                    KPP_taveFreq, and freq in pkg/diagnostics are converted to exact
C                    calendar months and years.  Requires pkg/cal.
      COMMON /CALENDAR_L/
     &               calendarDumps,
     &               usingModelCalendar,
     &               usingNoLeapYearCal,
     &               usingJulianCalendar,
     &               usingGregorianCalendar
      LOGICAL calendarDumps
      LOGICAL usingModelCalendar
      LOGICAL usingNoLeapYearCal
      LOGICAL usingJulianCalendar
      LOGICAL usingGregorianCalendar

C     theCalendar :: type of calendar to use; available:
C                    'model', 'gregorian' or 'noLeapYear'.
C     dayOfWeek   :: Week day number one is the week day of refDate.
C                    For the Gregorian calendar this is Friday, 15-Oct-1582.
C     monthOfYear :: Both available calendars are assumed to have twelve
C                    months.
      COMMON /CALENDAR_C/
     &                     theCalendar,
     &                     dayOfWeek,
     &                     monthOfYear
      CHARACTER*(20) theCalendar
      CHARACTER*(3) dayOfWeek(7)
      CHARACTER*(3) monthOfYear(nMonthYear)


c     == routine arguments ==

      integer iday
      integer firststep
      integer laststep
      integer nsteps
      integer mythid

c     == local variables ==

      integer ierr
      integer mdstep
      integer numdays
      integer numsteps
      integer frac1
      integer frac2
      integer frac3
      integer frac4
      integer fullsteps
      integer firstyear
      integer firstmonth
      integer firstday
      integer lyfirst
      integer startsecs
      integer lastday
      integer endsecs

c     == external ==

      external cal_IntDays
      integer  cal_IntDays
      external cal_IsLeap
      integer  cal_IsLeap

c     == end of interface ==

      numdays    = cal_IntDays( mythid )
      lyfirst    = cal_IsLeap( firstyear, mythid )

      mdstep     = int(modelstep)

      firstyear  = modelstartdate(1)/10000
      firstmonth = mod(modelstartdate(1)/100,100)
      firstday   = mod(modelstartdate(1),100)
      lastday    = mod(modelenddate(1),100)

      startsecs  = (modelstartdate(2)/10000)*secondsperhour +
     &             mod(modelstartdate(2)/100,100)*secondsperminute +
     &             mod(modelstartdate(2),100)
      endsecs    = (modelenddate(2)/10000)*secondsperhour +
     &             mod(modelenddate(2)/100,100)*secondsperminute +
     &             mod(modelenddate(2),100)

      if ( numdays .eq. 1 ) then
        if ( iday .eq. firstday ) then
c--       Get the number of steps in the first day.
          if ( firstday .eq. lastday ) then
            firststep = 1
            laststep  = modelintsteps
          else if ( mod(firstday+1,ndaymonth(firstmonth,lyfirst)) .eq.
     &              lastday ) then
c--         This can only happen if we end at midnight of the next day.
            if ( modelenddate(2) .eq. 0 ) then
              firststep = 1
              laststep  = modelintsteps
c--           Note: This holds only if modelenddate was determined
c--                 such that it coincides with the model final time.
            else
c--           We do not end at midnight of the first day of
c--           the next month.
              ierr = 2604
              call cal_PrintError( ierr, mythid )
              stop ' stopped in cal_StepsForDay.'
            endif
          else
c--         The first and the last day are inconsistent with iday.
            ierr = 2603
            call cal_PrintError( ierr, mythid )
            stop ' stopped in cal_StepsForDay.'
          endif
        else
c--       The variables numdays and iday are inconsistent;
c--       ( iday .gt. numdays ).
          ierr = 2602
          call cal_PrintError( ierr, mythid )
          stop ' stopped in cal_StepsForDay.'
        endif

      else if ( numdays .gt. 1 ) then
c--     More than one day of integration.
        if ( iday .eq. 1 ) then
          firststep = 1
          laststep  = int((secondsperday - startsecs)/mdstep)
        else if ( ( iday .gt. 1      ) .and.
     &            ( iday .lt. numdays) ) then
c--       Somewhere between first and last month.
c--       The first steps in iday.
          fullsteps = int((secondsperday - startsecs)/mdstep)
          numsteps  = fullsteps
c--       What is left in the first day (frac1).
          frac1     = (secondsperday - startsecs) - fullsteps*mdstep
          fullsteps = int(secondsperday/modelstep)
c--       What is left in a complete day (frac2).
          frac2     = secondsperday - fullsteps*mdstep
c--       What is left up to the current day (frac3).
          frac3     = frac1 + frac2*(iday - 1)
          numsteps  = numsteps + (iday - 1)*fullsteps +
     &                frac3/mdstep
          laststep  = numsteps
          firststep = laststep - secondsperday/mdstep + 1

        else if ( iday .eq. numdays ) then
c--       The last day of integration.
c--       The first step in iday.
          fullsteps = int((secondsperday - startsecs)/mdstep)
          numsteps  = fullsteps
c--       What is left in the first day (frac1).
          frac1     = (secondsperday - startsecs) - fullsteps*mdstep
          fullsteps = int(secondsperday/modelstep)
c--       What is left in a complete day (frac2).
          frac2     = secondsperday - fullsteps*mdstep
c--       What is left up to the day before the last (frac3).
          frac3     = frac1 + frac2*(iday - 2)
          numsteps  = numsteps + (iday - 2)*fullsteps
c--       The last step in iday.
          if ( modelenddate(2) .eq. 0 ) then
c--         This can only happen if we end at midnight of the next day.
            laststep  = numsteps + fullsteps
            firststep = numsteps + 1
c--         Note: There should be no fraction left
c--               ( mod(frac3,mdstep) = frac3/mdstep ) if modelenddate
c--               is based on an integral number of timesteps.
          else
            frac4     = frac3 + endsecs
            numsteps  = numsteps + frac4/mdstep
            laststep  = numsteps

c--         Note: There should be no fraction left
c--               ( mod(frac4,mdstep = frac4/mdstep ) if modelenddate
c--               is based on an integral number of timesteps.

            firststep = laststep - endsecs/mdstep + 1
          endif
        else
c--       The variables iday and numdays are inconsistent.
          ierr = 2605
          call cal_PrintError( ierr, mythid )
          stop ' stopped in cal_DaysForMonth.'
        endif
      else
c--     The number of days to integrate is wrong; check cal_IntDays.
        ierr = 2601
        call cal_PrintError( ierr, mythid )
        stop ' stopped in cal_StepsForDay.'
      endif

c--   The number of days to integrate in the given month.
      nsteps = laststep - firststep + 1

      return
      end

