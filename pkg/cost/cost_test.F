#include "COST_OPTIONS.h"

      subroutine cost_test( myThid )
C     *==========================================================*
C     | subroutine cost_test
C     | o this routine computes the cost function for the tiles
C     |   of this processor
C     *==========================================================*
C     | Notes
C     *==========================================================*
      IMPLICIT NONE

C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"

#include "cost.h"

C     == Routine arguments ==
C     myThid - Thread number for this instance of the routine.
      integer myThid

#ifdef ALLOW_COST_TEST
C     == Local variables
      integer bi, bj
      _RL thetaRef
      integer i, j, k

      LOGICAL exst
      _RL tmp, thetaTmp, volwght, locfc
      _RL dummyRL(1)
      _RS mask_ws        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL volWghtFld     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      INTEGER ILNBLNK, il

      if ( useCubedSphereExchange) then
         iLocOut =  7
         jLocOut =  28
         kLocOut =  1
      else
         iLocOut =  80
         jLocOut =  30
         kLocOut =  1
      endif

ce    some reference temperature
c      thetaRef = 24.0 _d 0
      thetaRef = 0.0 _d 0

C--   Calculate mask for tracer cells  (0 => land, 1 => water)
      k=1

C--   Read file
      exst=.FALSE.
      il = ILNBLNK(cost_mask_file)
      INQUIRE(file = cost_mask_file(1:il), exist = exst)
#ifdef ALLOW_MDSIO
      IF ( exst ) THEN
       CALL MDS_READ_FIELD(
     I      cost_mask_file(1:il), readBinaryPrec, .FALSE., 'RS',
     I      Nr, 1, Nr,
     O      dummyRL, mask_ws,
     I      1, myThid )
      ELSE
       WRITE(standardmessageunit,*)
     &      'S/R COST_TEST: ', cost_mask_file(1:il), ' does not exist'
      ENDIF
#else
      cost_mean_salt.F: Oops!
#endif

C  first compute volume weights
      tmp  = 0. _d 0
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        volWghtFld(:,:,:,bi,bj) = 0. _d 0
        DO k=1,Nr
         DO j=1,sNy
          DO i=1,sNx
           volwght = mask_ws(i,j,k,bi,bj)
     &          *drF(k)*hFacC(i,j,k,bi,bj)*rA(i,j,bi,bj)
CML           volwght = 1. _d 0
           volWghtFld(i,j,k,bi,bj) = volwght
           tmp = tmp + volwght
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      _GLOBAL_SUM_RL( tmp , myThid )
      IF ( tmp.GT.0. ) tmp = 1. _d 0 / tmp
C
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        locfc = 0. _d 0
        DO k=1,Nr
         DO j=1,sNy
          DO i=1,sNx
           thetaTmp = cMeanTheta(i,j,k,bi,bj) - thetaRef
CML           thetaTmp = theta(i,j,k,bi,bj) - thetaRef
C     scale with tmp to reduce truncation error?
           locfc = locfc + volWghtFld(i,j,k,bi,bj) * thetaTmp
CML           print*, 'mean theta: ',cMeanTheta(i,j,k,bi,bj),
CML     &          theta(i,j,k,bi,bj), thetaTmp
CML           print*, 'weddell sea mask:', mask_ws(i,j,k,bi,bj)
CML           print*, 'weighted volume', volwghtfld(i,j,k,bi,bj)
CML           print*, 'mean of tile:',locfc
          ENDDO
         ENDDO
        ENDDO
        objf_test(bi,bj) = locfc * tmp
       ENDDO
      ENDDO

C     diagnostics
      CALL WRITE_FLD_XYZ_RL('maskWSVol',' ',volWghtFld,-1,myThid)
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        write(standardmessageunit,*)
     &       'ml-objf_test     =',objf_test(bi,bj),
     &       tmp,startTime,endTime,lastInterval
       ENDDO
      ENDDO

#endif  /* ALLOW_COST_TEST */

      RETURN
      END
