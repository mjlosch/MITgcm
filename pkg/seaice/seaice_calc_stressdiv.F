#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C--  File seaice_calc_stressdiv.F
C--   Contents
C--   o SEAICE_CALC_STRESSDIV
C--   o SEAICE_CALC_STRESS

CBOP
C !ROUTINE: SEAICE_CALC_STRESSDIV
C !INTERFACE: ==========================================================
      SUBROUTINE SEAICE_CALC_STRESSDIV(
     I     e11l, e22l, e12l, pressl, zetal, etal, etaZl,
     O     stressDivX, stressDivY,
     I     bi, bj, myTime, myIter, myThid )

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE SEAICE_CALC_STRESSDIV
C     | o compute ice internal stress divergence
C     |
C     | Martin Losch, May 2017, Martin.Losch@awi.de
C     *===========================================================*
C \ev

C !USES: ===============================================================
      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     e11/e22/e12 :: strain rate tensor components
C     press  :: maximal compressive strength
C     zeta   :: bulk viscosity
C     eta    :: shear viscosity
C     etaZ   :: shear viscosity at vorticity points
C     stressDivX/Y :: x/y-component of stress divergence
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
      _RL e11l  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22l  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12l  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pressl(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zetal (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etal  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etaZl (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef SEAICE_CGRID
C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j       :: inner loop counters
      _RL sig11(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig22(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER i, j
      _RL hFacU, hFacV
      _RL e12W, e12E, e12S, e12N

      CALL SEAICE_CALC_STRESS(
     I     e11l, e22l, e12l, pressl, zetal, etal, etaZl,
     O     sig11, sig22, sig12,
     I     bi, bj, myTime, myIter, myThid )

C     compute divergence of stress tensor

      DO J=1,sNy
       DO I=1,sNx
        stressDivX(i,j,bi,bj) =
     &       ( sig11(i  ,j  ) * _dyF(i  ,j,bi,bj)
     &       - sig11(i-1,j  ) * _dyF(i-1,j,bi,bj)
     &       + sig12(i  ,j+1) * _dxV(i,j+1,bi,bj)
     &       - sig12(i  ,j  ) * _dxV(i,j  ,bi,bj)
     &       ) * recip_rAw(i,j,bi,bj)
        stressDivY(i,j,bi,bj) =
     &       ( sig22(i  ,j  ) * _dxF(i,j  ,bi,bj)
     &       - sig22(i  ,j-1) * _dxF(i,j-1,bi,bj)
     &       + sig12(i+1,j  ) * _dyU(i+1,j,bi,bj)
     &       - sig12(i  ,j  ) * _dyU(i  ,j,bi,bj)
     &       ) * recip_rAs(i,j,bi,bj)
       ENDDO
      ENDDO
      IF ( SEAICE_no_slip ) THEN
       DO J=1,sNy
        DO I=1,sNx
         hFacU = SIMaskU(i,j,bi,bj) - SIMaskU(i,j-1,bi,bj)
         e12E = (1. - SImaskU(i+1,j,bi,bj)*SImaskU(i+1,j-1,bi,bj))
     &        * (- vIce(i,j,bi,bj) *  _recip_dxV(i+1,j,bi,bj))
     &        + ( uIce(i,j,bi,bj) + uIce(i,j-1,bi,bj) )
     &        * _recip_dyU(i,j,bi,bj) * hFacU
         hFacU = SIMaskU(i,j+1,bi,bj) - SIMaskU(i,j,bi,bj)
         e12W = (1. - SImaskU(i  ,j,bi,bj)*SImaskU(i  ,j-1,bi,bj))
     &        *    vIce(i,j,bi,bj) *  _recip_dxV(i  ,j,bi,bj)
     &        + ( uIce(i,j,bi,bj) + uIce(i,j+1,bi,bj) )
     &        * _recip_dyU(i,j+1,bi,bj) * hFacU
         hFacV = SIMaskV(i,j,bi,bj) - SIMaskV(i-1,j,bi,bj)
         e12S = (1. - SImaskV(i,j  ,bi,bj)*SImaskV(i-1,j  ,bi,bj))
     &        *    uIce(i,j,bi,bj) *  _recip_dyU(i,j  ,bi,bj)
     &        + ( vIce(i,j,bi,bj) + vIce(i-1,j,bi,bj) )
     &        * _recip_dxV(i,j,bi,bj) * hFacV
         hFacV = SIMaskV(i+1,j,bi,bj) - SIMaskV(i,j,bi,bj)
         e12N = (1. - SImaskV(i,j+1,bi,bj)*SImaskV(i-1,j+1,bi,bj))
     &        * (- uIce(i,j,bi,bj) *  _recip_dyU(i,j+1,bi,bj))
     &        + ( vIce(i,j,bi,bj) + vIce(i+1,j,bi,bj) )
     &        * _recip_dxV(i,j+1,bi,bj) * hFacV
CML       if (myiter.eq.1) print *, 'ml-e12: ', i,j,e12S, e12N, e12W, e12E
         stressDivX(i,j,bi,bj) = stressDivX(i,j,bi,bj) +
     &        ( 2.*etaZl(i  ,j+1,bi,bj) * e12N * _dxV(i,j+1,bi,bj)
     &        - 2.*etaZl(i  ,j  ,bi,bj) * e12S * _dxV(i,j  ,bi,bj)
     &        ) * recip_rAw(i,j,bi,bj)
         stressDivY(i,j,bi,bj) = stressDivY(i,j,bi,bj) +
     &        ( 2.*etaZl(i+1,j  ,bi,bj) * e12E * _dyU(i+1,j,bi,bj)
     &        - 2.*etaZl(i  ,j  ,bi,bj) * e12W * _dyU(i  ,j,bi,bj)
     &        ) * recip_rAs(i,j,bi,bj)
        ENDDO
       ENDDO
      ENDIF

      RETURN
      END

CBOP
C !ROUTINE: SEAICE_CALC_STRESS
C !INTERFACE: ==========================================================
      SUBROUTINE SEAICE_CALC_STRESS(
     I     e11, e22, e12, press, zeta, eta, etaZ,
     O     sig11, sig22, sig12,
     I     bi, bj, myTime, myIter, myThid )

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE SEAICE_CALC_STRESS
C     | o compute ice internal stress according to rheology
C     |
C     | Martin Losch, May 2017, Martin.Losch@awi.de
C     *===========================================================*
C \ev

C !USES: ===============================================================
      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
C     e11/e22/e12 :: strain rate tensor components
C     press  :: maximal compressive strength
C     zeta   :: bulk viscosity
C     eta    :: shear viscosity
C     etaZ   :: shear viscosity at vorticity points
C     sig11/sig22/sig12   :: stress tensor components
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER bi, bj
      _RL e11  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zeta (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL eta  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etaZ (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sig11(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig22(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CEOP

C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j           :: inner loop counters
C     eplus, eminus :: convenient abbreviations for e11+e22, e11-e22
      INTEGER i, j
      _RL     eplus, eminus

C     compute components of stress tensor from current strain rate tensor
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        sig11(i,j) = 0. _d 0
        sig22(i,j) = 0. _d 0
        sig12(i,j) = 0. _d 0
       ENDDO
      ENDDO

      DO j=0,sNy
       DO i=0,sNx
        eplus = e11(i,j,bi,bj) + e22(i,j,bi,bj)
        eminus= e11(i,j,bi,bj) - e22(i,j,bi,bj)
        sig11(i,j) = zeta(i,j,bi,bj)*eplus + eta(i,j,bi,bj)*eminus
     &       - 0.5 _d 0 * press(i,j,bi,bj)
        sig22(i,j) = zeta(i,j,bi,bj)*eplus - eta(i,j,bi,bj)*eminus
     &       - 0.5 _d 0 * press(i,j,bi,bj)
       ENDDO
      ENDDO

      DO j=1,sNy+1
       DO i=1,sNx+1
        sig12(i,j) = 2. _d 0 * e12(i,j,bi,bj) * etaZ(i,j,bi,bj)
       ENDDO
      ENDDO
#endif /* SEAICE_CGRID */
      RETURN
      END
