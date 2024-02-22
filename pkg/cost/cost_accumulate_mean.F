#include "COST_OPTIONS.h"

      subroutine cost_accumulate_mean( myThid )
C     *==========================================================*
C     | subroutine cost_accumulate_mean                          |
C     | o accumulate mean state for cost evalualtion             |
C     *==========================================================*
C     |                                                          |
C     *==========================================================*
      IMPLICIT NONE

C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "FFIELDS.h"

#include "cost.h"

C     == Routine arguments ==
C     myThid - Thread number for this instance of the routine.
      integer bi, bj
      integer myThid

#ifdef ALLOW_COST
C     == Local variables
      _RL thetaRef
      _RL dTfac

      integer i, j, k
      integer ig, jg
      integer itlo,ithi
      integer jtlo,jthi

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

      dTfac = deltaTClock / lastInterval

C--   Calculate cost function on tile of this instance
      do bj = jtlo,jthi
        do bi = itlo,ithi
          do k = 1, Nr
            do j=1,sNy
              do i=1,sNx
                cMeanTheta(i,j,k,bi,bj) = cMeanTheta(i,j,k,bi,bj)
     &                + theta(i,j,k,bi,bj) * dTfac
                cMeanSalt(i,j,k,bi,bj) = cMeanSalt(i,j,k,bi,bj)
     &                + salt(i,j,k,bi,bj)  * dTfac
                cMeanUVel(i,j,k,bi,bj) = cMeanUVel(i,j,k,bi,bj)
     &               + uVel(i,j,k,bi,bj)   * dTfac
                cMeanVVel(i,j,k,bi,bj) = cMeanVVel(i,j,k,bi,bj)
     &               + vVel(i,j,k,bi,bj)   * dTfac

                cMeanThetaUVel(i,j,k,bi,bj) =
     &               cMeanThetaUVel(i,j,k,bi,bj)
     &               + (theta(i,j,k,bi,bj)+theta(i-1,j,k,bi,bj))
     &                 *0.5 _d 0*uvel(i,j,k,bi,bj)
     &                 *maskW(i,j,k,bi,bj)*maskC(i,j,k,bi,bj)
     &                 *dTfac
                cMeanThetaVVel(i,j,k,bi,bj) =
     &               cMeanThetaVVel(i,j,k,bi,bj)
     &               + (theta(i,j,k,bi,bj)+theta(i,j-1,k,bi,bj))
     &                 *0.5 _d 0*vvel(i,j,k,bi,bj)
     &                 *maskS(i,j,k,bi,bj)*maskC(i,j,k,bi,bj)
     &                 *dTfac

                cMeanSaltUVel(i,j,k,bi,bj) =
     &               cMeanSaltUVel(i,j,k,bi,bj)
     &               + (salt(i,j,k,bi,bj)+salt(i-1,j,k,bi,bj))
     &                 *0.5 _d 0*uvel(i,j,k,bi,bj)
     &                 *maskW(i,j,k,bi,bj)*maskC(i,j,k,bi,bj)
     &                 *dTfac
                cMeanSaltVVel(i,j,k,bi,bj) =
     &               cMeanSaltVVel(i,j,k,bi,bj)
     &               + (salt(i,j,k,bi,bj)+salt(i,j-1,k,bi,bj))
     &                 *0.5 _d 0*vvel(i,j,k,bi,bj)
     &                 *maskS(i,j,k,bi,bj)*maskC(i,j,k,bi,bj)
     &                 *dTfac
              end do
            end do
          end do
          do j=1,sNy
            do i=1,sNx
              cMeanEtaTot(i,j,bi,bj) = cMeanEtaTot(i,j,bi,bj)
     &            + ( etaN(i,j,bi,bj)
     &              + sIceLoad(i,j,bi,bj) * recip_rhoConst
     &              ) * dTfac
            end do
          end do
        end do
      end do

#endif

      RETURN
      END
