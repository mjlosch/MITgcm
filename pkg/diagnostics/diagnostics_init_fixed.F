#include "DIAG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C     !ROUTINE: DIAGNOSTICS_INIT_FIXED

C     !INTERFACE:
      SUBROUTINE DIAGNOSTICS_INIT_FIXED(myThid)

C     !DESCRIPTION:
C     finish setting up the list of available diagnostics and
C     prepare for storing selected diagnostics and statistics-diags.

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "DIAGNOSTICS_SIZE.h"
#include "DIAGNOSTICS.h"
#ifdef ALLOW_FIZHI
#include "PARAMS.h"
#endif

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
c#ifdef ALLOW_DIVIDED_ADJOINT
      LOGICAL dummyBoolean
c#endif

C--   Set number of levels for all available diagnostics
C     (cannot add diags to list anymore after this call)
      CALL DIAGNOSTICS_SET_LEVELS( myThid )

C--   Calculate pointers for diagnostics set to non-zero frequency
      CALL DIAGNOSTICS_SET_POINTERS( myThid )

C--   Set parameters and variables used in post-processing diagnostics
      CALL DIAGNOSTICS_SET_CALC( myThid )

C--   Define region-mask for regional statistics diagnostics
      CALL DIAGSTATS_SET_REGIONS( myThid )

C--   Calculate pointers for statistics-diags set to non-zero frequency
      CALL DIAGSTATS_SET_POINTERS( myThid )

#ifdef ALLOW_FIZHI
      IF ( useFIZHI ) THEN
        CALL FIZHI_DIAGALARMS( myThid )
      ENDIF
#endif
c#ifdef ALLOW_DIVIDED_ADJOINT
      CALL DIAGS_TRACK_DIVA( dummyBoolean, -2, myThid )
c#endif

      RETURN
      END
