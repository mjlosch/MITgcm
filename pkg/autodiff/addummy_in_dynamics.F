#include "AUTODIFF_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#include "AD_CONFIG.h"

CBOP
C     !ROUTINE: ADDUMMY_IN_DYNAMICS
C     !INTERFACE:
      SUBROUTINE ADDUMMY_IN_DYNAMICS( myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE ADDUMMY_IN_DYNAMICS                           |
C     *==========================================================*
C     Extract adjoint variable from TAMC/TAF-generated
C     adjoint common blocks, contained in adcommon.h
C     and write fields to file;
C     Make sure common blocks in adcommon.h are up-to-date
C     w.r.t. current adjoint code.
C     *==========================================================*
C     | SUBROUTINE ADDUMMY_IN_DYNAMICS                           |
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "AUTODIFF_PARAMS.h"
#ifdef ALLOW_MNC
#include "MNC_PARAMS.h"
#endif
c#include "GRID.h"
#ifdef ALLOW_AUTODIFF_MONITOR
# include "adcommon.h"
# include "AUTODIFF.h"
# ifdef ALLOW_DIAGNOSTICS
#  include "DIAGNOSTICS_P2SHARE.h"
# endif
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myTime :: time counter for this thread
C     myIter :: iteration counter for this thread
C     myThid :: Thread number for this instance of the routine.
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

#if (defined (ALLOW_ADJOINT_RUN) || defined (ALLOW_ADMTLM))
# ifdef ALLOW_AUTODIFF_MONITOR

C     !FUNCTIONS:
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C     !LOCAL VARIABLES:
c     == local variables ==
C     suff   :: Hold suffix part of a filename
C     msgBuf :: Error message buffer
      CHARACTER*(10) suff
c     CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RS dumRS(1)
      LOGICAL doExch
      LOGICAL doDump
CEOP

C--- Do dump?
      doDump = DIFFERENT_MULTIPLE( adjDumpFreq, myTime, deltaTClock )

C--- Test for adexch based on dumpAdVarExch flag
      doExch = doDump
#  ifdef ALLOW_DIAGNOSTICS
      doExch = doDump .OR. useDiag4AdjOutp
#  endif
      doExch = doExch .AND. ( dumpAdVarExch.EQ.1 )

C--- If necessary, do the exchanges
      IF ( doExch ) THEN
#  ifdef ALLOW_OPENAD
C--   need to all the correct OpenAD EXCH S/R ; left empty for now
#  else /* ALLOW_OPENAD */

#   ifdef AUTODIFF_TAMC_COMPATIBILITY
C
#   else /* ndfef AUTODIFF_TAMC_COMPATIBILITY */
#    ifdef ALLOW_AUTODIFF_MONITOR_DIAG
       CALL ADEXCH_3D_RL( adRhoInSitu, Nr, myThid )
       CALL ADEXCH_3D_RL( adTotPhiHyd, Nr, myThid )
#    endif
#   endif /* AUTODIFF_TAMC_COMPATIBILITY */
#  endif
      ENDIF


      IF ( doDump ) THEN
C     Increment ad dump record number (used only if dumpAdByRec is true)
        dumpAdRecMn = dumpAdRecMn + 1
        CALL TIMER_START('I/O (WRITE)        [ADJOINT LOOP]', myThid )
      ENDIF
#  ifndef ALLOW_OPENAD
C========================
C--  TAF AD-variables:
C========================
      CALL DUMP_ADJ_XYZ(dumRS, adRhoInSitu, 'ADJrho  ', 'ADJrhoInSitu.',
     &                  12, doDump, dumpAdRecMn, myTime, myIter,myThid)
      CALL DUMP_ADJ_XYZ(dumRS, adTotPhiHyd, 'ADJtPhiH', 'ADtotPhiHyd.',
     &                  12, doDump, dumpAdRecMn, myTime, myIter,myThid)
#  else
      CALL DUMP_ADJ_XYZ(dumRS, rhoInSitu%d, 'ADJrho  ', 'ADJrhoInSitu.',
     &                  12, doDump, dumpAdRecMn, myTime, myIter,myThid)
      CALL DUMP_ADJ_XYZ(dumRS, totPhiHyd%d, 'ADJtPhiH', 'ADtotPhiHyd.',
     &                  12, doDump, dumpAdRecMn, myTime, myIter,myThid)
#  endif

      IF ( doDump ) THEN
       CALL TIMER_STOP( 'I/O (WRITE)        [ADJOINT LOOP]', myThid )
      ENDIF

# endif /* ALLOW_AUTODIFF_MONITOR */
#endif /* ALLOW_ADJOINT_RUN */

      RETURN
      END
