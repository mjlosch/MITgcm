#include "COST_OPTIONS.h"
c#ifdef ALLOW_CTRL
c# include "CTRL_OPTIONS.h"
c#endif

      SUBROUTINE COST_MEAN_SALT( myThid )
C     *==========================================================*
C     | SUBROUTINE COST_TEMP
C     | o the subroutine computes the sum of the squared errors
C     |   relatively to the Levitus climatology
C     *==========================================================*

       IMPLICIT NONE

C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "cost.h"
#include "ctrl_weights.h"

C     ======== Routine arguments ======================
C     myThid - Thread number for this instance of the routine.
      INTEGER myThid

#ifdef ALLOW_COST_MEAN_SALT
C     ========= Local variables =========================
      INTEGER i, j, k
      INTEGER bi, bj
      INTEGER Nk
      _RL locfc,tmp,volwght,saltTmp
      _RS mask_arc_oc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,   nSx,nSy)
      _RL volWghtFld (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer itlo,ithi
      integer jtlo,jthi

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

      Nk = 2
C Read annual mean Levitus temperature

      CALL READ_FLD_XY_RS('mask_arctic_ocean.data',
     &                    ' ',mask_arc_oc,0,myThid)

CML      CALL WRITE_FLD_XY_RS('maskArcticTest',' ',mask_arc_oc,-1,myThid)

C  first compute volume weights
      tmp  = 0. _d 0
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        volWghtFld(:,:,:,bi,bj) = 0. _d 0
        DO k=1,Nr
         DO j=1,sNy
          DO i=1,sNx
           saltTmp = cMeanSalt(i,j,k,bi,bj)
CML           saltTmp = salt(i,j,k,bi,bj)
           IF ( saltTmp .LE. 34. _d 0 ) THEN
            volwght = mask_arc_oc(i,j,bi,bj)
     &          *drF(k)*hFacC(i,j,k,bi,bj)*rA(i,j,bi,bj)
CML     &          *(1.-1./(1.+exp(-20. _d 0 *(saltTmp-35. _d 0)))
           ELSE
            volwght = 0. _d 0
           ENDIF
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
           saltTmp = cMeanSalt(i,j,k,bi,bj) - 35. _d 0
CML           saltTmp = salt(i,j,k,bi,bj) - 35. _d 0
C     scale with tmp to reduce truncation error?
           locfc = locfc + volWghtFld(i,j,k,bi,bj) * saltTmp ! * tmp
CML           print*, 'mean salinity: ',cMeanSalt(i,j,k,bi,bj)
CML           print*, 'arctic coean mask:', mask_arc_oc(i,j,bi,bj)
CML           print*, 'logistic equation:'
CML           print*, (1-1/(1+exp(-20*(cMeanSalt(i,j,k,bi,bj)-35))))
CML           print*, 'weighted volume', volwght
CML           print*, 'mean of tile:',locfc
          ENDDO
         ENDDO
        ENDDO
        objf_mean_salt(bi,bj) = locfc
       ENDDO
      ENDDO

C     diagnostics
      CALL WRITE_FLD_XYZ_RL('maskArcticVol',' ',volWghtFld,-1,myThid)
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        write(standardmessageunit,*)
     &       'ml-objf_mean_salt =',objf_mean_salt(bi,bj),
     &       tmp,startTime,endTime
       ENDDO
      ENDDO

#endif /* ALLOW_COST_MEAN_SALT */
      RETURN
      END
