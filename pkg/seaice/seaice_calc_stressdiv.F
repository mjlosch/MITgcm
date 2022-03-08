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
     I     e11, e22, e12, press, zeta, eta, etaZ,
     O     stressDivergenceX, stressDivergenceY,
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

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     e11/e22/e12 :: strain rate tensor components
C     press  :: maximal compressive strength
C     zeta   :: bulk viscosity
C     eta    :: shear viscosity
C     etaZ   :: shear viscosity at vorticity points
C     stressDivergenceX/Y :: x/y-component of stress divergence
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
      _RL e11  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zeta (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL eta  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etaZ (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
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

      print *, 'ml-hallo', myIter
      IF ( SEAICEuseVelStress .AND. myIter .GE. 0 ) THEN
C     special treatment of stress in this case:
C     e11=uIce, e22=vIce, e12=duIce, press=dvIce, zeta=press,
C     eta=press0, etaZ=deltaC
       CALL SEAICE_CALC_STRESS_VELSTRESS(
     I      e11, e22, e12, press, zeta, eta, etaZ,
     O      sig11, sig22, sig12,
     I      bi, bj, myTime, myIter, myThid )
      ELSE
       CALL SEAICE_CALC_STRESS(
     I      e11, e22, e12, press, zeta, eta, etaZ,
     O      sig11, sig22, sig12,
     I      bi, bj, myTime, myIter, myThid )
      ENDIF

C     compute divergence of stress tensor

      DO j=1,sNy
       DO i=1,sNx
        stressDivergenceX(i,j,bi,bj) =
     &       ( sig11(i  ,j  ) * _dyF(i  ,j,bi,bj)
     &       - sig11(i-1,j  ) * _dyF(i-1,j,bi,bj)
     &       + sig12(i  ,j+1) * _dxV(i,j+1,bi,bj)
     &       - sig12(i  ,j  ) * _dxV(i,j  ,bi,bj)
     &       ) * recip_rAw(i,j,bi,bj)
        stressDivergenceY(i,j,bi,bj) =
     &       ( sig22(i  ,j  ) * _dxF(i,j  ,bi,bj)
     &       - sig22(i  ,j-1) * _dxF(i,j-1,bi,bj)
     &       + sig12(i+1,j  ) * _dyU(i+1,j,bi,bj)
     &       - sig12(i  ,j  ) * _dyU(i  ,j,bi,bj)
     &       ) * recip_rAs(i,j,bi,bj)
       ENDDO
      ENDDO

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

      RETURN
      END

CBOP
C !ROUTINE: SEAICE_CALC_STRESS_VELSTRESS
C !INTERFACE: ==========================================================
      SUBROUTINE SEAICE_CALC_STRESS_VELSTRESS(
     I     uLoc, vLoc, duIce, dvIce, pr, pr0, delC,
     O     sig11, sig22, sig12,
     I     bi, bj, myTime, myIter, myThid )
C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE SEAICE_CALC_STRESS_VELSTRESS
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
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
C     u/vLoc :: ice velcoities
C     du/vLoc:: ice velcoity perturbations
C     sig11/sig22/sig12   :: stress tensor components
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER bi, bj
      _RL uLoc (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vLoc (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL duIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dvIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pr   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pr0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL delC (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sig11(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig22(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CEOP

C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j           :: inner loop counters
C     eplus, eminus :: convenient abbreviations for e11+e22, e11-e22
      INTEGER i, j
      _RL recip_ecc
      _RL pp, pm, dpp, dpm
      _RL pressZ, sumNorm

      recip_ecc = 1./SEAICE_eccen

      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        sig11(i,j) = 0. _d 0
        sig22(i,j) = 0. _d 0
        sig12(i,j) = 0. _d 0
CML        e11p(i,j,bi,bj) = 0.
CML        e22p(i,j,bi,bj) = 0.
CML        e12p(i,j,bi,bj) = 0.
       ENDDO
      ENDDO

      CALL SEAICE_CALC_DPI(
CML     I     e11, e22, e12, delC,
CML     O     dpi11, dpi22, dpi12,
     I     bi, bj, myTime, myIter, myThid )
C
      DO j=0,sNy
       DO i=0,sNx
        pp  = ( pi11(i,j,bi,bj)+ pi22(i,j,bi,bj))
        pm  = ( pi11(i,j,bi,bj)- pi22(i,j,bi,bj))
        dpp = (dpi11(i,j,bi,bj)+dpi22(i,j,bi,bj))
        dpm = (dpi11(i,j,bi,bj)-dpi22(i,j,bi,bj))
        sig11(i,j) = (              0.5 _d 0 *(  pp + recip_ecc * pm )
     &       + SEAICE_JFNKepsilon * 0.5 _d 0 *( dpp + recip_ecc *dpm )
     &       ) * pr0(i,j,bi,bj)*(ONE+tensileStrFac(i,j,bi,bj))
        sig22(i,j) = (              0.5 _d 0 *(  pp - recip_ecc * pm )
     &       + SEAICE_JFNKepsilon * 0.5 _d 0 *( dpp - recip_ecc *dpm )
     &       ) * pr0(i,j,bi,bj)*(ONE+tensileStrFac(i,j,bi,bj))
        sig11(i,j) = sig11(i,j) - 0.5 _d 0 * pr(i,j,bi,bj)
        sig22(i,j) = sig22(i,j) - 0.5 _d 0 * pr(i,j,bi,bj)
       ENDDO
      ENDDO
      DO j=1,sNy+1
       DO i=1,sNx+1
        sumNorm = HEFFM(i,j,  bi,bj)+HEFFM(i-1,j,  bi,bj)
     &          + HEFFM(i,j-1,bi,bj)+HEFFM(i-1,j-1,bi,bj)
        IF ( sumNorm.GT.0. _d 0 ) sumNorm = 1. _d 0 / sumNorm
        pressZ = sumNorm *
     &       ( pr0(i,  j,  bi,bj)*(ONE+tensileStrFac(i,  j,  bi,bj))
     &       + pr0(i-1,j-1,bi,bj)*(ONE+tensileStrFac(i-1,j-1,bi,bj))
     &       + pr0(i-1,j,  bi,bj)*(ONE+tensileStrFac(i-1,j,  bi,bj))
     &       + pr0(i,  j-1,bi,bj)*(ONE+tensileStrFac(i,  j-1,bi,bj))
     &       )
        sig12(i,j) = recip_ecc * ( pi12(i,j,bi,bj)
     &       + SEAICE_JFNKepsilon * dpi12(i,j,bi,bj)
     &       ) * pressZ
       ENDDO
      ENDDO
      RETURN
      END

CBOP
C !ROUTINE: SEAICE_CALC_DPI
C !INTERFACE: ==========================================================
      SUBROUTINE SEAICE_CALC_DPI(
CML     I     e110, e220, e120, delC,
CML     O     ddpi11, ddpi22, ddpi12,
     I     bi, bj, myTime, myIter, myThid )
C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE SEAICE_CALC_DPI
C     | o compute tensor dpi
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
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
C     u/vLoc :: ice velcoities
C     du/vLoc:: ice velcoity perturbations
C     sig11/sig22/sig12   :: stress tensor components
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER bi, bj
CML      _RL delC (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL e110(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL e220(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL e120(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL ddpi11(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL ddpi22(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL ddpi12(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CEOP

C !LOCAL VARIABLES: ====================================================
C     === Local variables ===
C     i,j           :: inner loop counters
C     eplus, eminus :: convenient abbreviations for e11+e22, e11-e22
      INTEGER i, j
      _RL recip_ecc, recip_esq, deltaReg, deltaRegSq, dsq
      _RL recip_delta, denom
      _RL ep, em, epp, emp
      _RL tau11, tau22, tau12
      _RL dtau11, dtau22, dtau12
      _RL pi12C, tau12C, dtau12C, c1Z, c2Z, sumNorm, dtp12C, dtt12C
      _RL c1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL c2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      recip_ecc = 1./SEAICE_eccen
      recip_esq = 1./SEAICE_eccen**2
      dsq = SEAICE_deltaMin**2

C     initialisze
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        dpi11(i,j,bi,bj) = 0.
        dpi22(i,j,bi,bj) = 0.
        dpi12(i,j,bi,bj) = 0.
       ENDDO
      ENDDO
C     first compute increment of tensor pi
      DO j=0,sNy+1
       DO i=0,sNx+1
#ifdef SEAICE_DELTA_SMOOTHREG
        deltaRegSq = deltaC(i,j,bi,bj)**2 + dsq
        deltaReg   = SQRT(deltaRegSq)
C        deltaReg = deltaC(i,j,bi,bj) + SEAICE_deltaMin
C        deltaRegSq = deltaReg**2
#else
        deltaReg   = MAX(deltaC(i,j,bi,bj), SEAICE_deltaMin)
        deltaRegSq = deltaReg**2
#endif
        recip_delta = 1. _d 0 / deltaReg
        ep     = e11 (i,j,bi,bj) + e22 (i,j,bi,bj)
        em     = e11 (i,j,bi,bj) - e22 (i,j,bi,bj)
        epp    = e11p(i,j,bi,bj) + e22p(i,j,bi,bj)
        emp    = e11p(i,j,bi,bj) - e22p(i,j,bi,bj)
        tau11  = 0.5*(ep  + recip_ecc*em)
        tau22  = 0.5*(ep  - recip_ecc*em)
        dtau11 = 0.5*(epp + recip_ecc*emp)
        dtau22 = 0.5*(epp - recip_ecc*emp)
        pi12C = SQRT( 0.25 _d 0 * recip_rA(i,j,bi,bj) *
     &       ( rAz(i  ,j  ,bi,bj)*pi12(i  ,j  ,bi,bj)**2
     &       + rAz(i+1,j  ,bi,bj)*pi12(i+1,j  ,bi,bj)**2
     &       + rAz(i  ,j+1,bi,bj)*pi12(i  ,j+1,bi,bj)**2
     &       + rAz(i+1,j+1,bi,bj)*pi12(i+1,j+1,bi,bj)**2 ) )
        tau12C = recip_ecc * SQRT( 0.25 _d 0 * recip_rA(i,j,bi,bj) *
     &       ( rAz(i  ,j  ,bi,bj)*e12(i  ,j  ,bi,bj)**2
     &       + rAz(i+1,j  ,bi,bj)*e12(i+1,j  ,bi,bj)**2
     &       + rAz(i  ,j+1,bi,bj)*e12(i  ,j+1,bi,bj)**2
     &       + rAz(i+1,j+1,bi,bj)*e12(i+1,j+1,bi,bj)**2 ) )
        dtau12C = recip_ecc * SQRT( 0.25 _d 0 * recip_rA(i,j,bi,bj) *
     &       ( rAz(i  ,j  ,bi,bj)*e12p(i  ,j  ,bi,bj)**2
     &       + rAz(i+1,j  ,bi,bj)*e12p(i+1,j  ,bi,bj)**2
     &       + rAz(i  ,j+1,bi,bj)*e12p(i  ,j+1,bi,bj)**2
     &       + rAz(i+1,j+1,bi,bj)*e12p(i+1,j+1,bi,bj)**2 ) )
        dtt12C = recip_esq * 0.25 _d 0 * recip_rA(i,j,bi,bj) *
     &   ( rAz(i  ,j  ,bi,bj)*e12p(i  ,j  ,bi,bj)*e12(i  ,j  ,bi,bj)
     &   + rAz(i+1,j  ,bi,bj)*e12p(i+1,j  ,bi,bj)*e12(i+1,j  ,bi,bj)
     &   + rAz(i  ,j+1,bi,bj)*e12p(i  ,j+1,bi,bj)*e12(i  ,j+1,bi,bj)
     &   + rAz(i+1,j+1,bi,bj)*e12p(i+1,j+1,bi,bj)*e12(i+1,j+1,bi,bj) )
        dtp12C = recip_ecc * 0.25 _d 0 * recip_rA(i,j,bi,bj) *
     &   ( rAz(i  ,j  ,bi,bj)*e12p(i  ,j  ,bi,bj)*pi12(i  ,j  ,bi,bj)
     &   + rAz(i+1,j  ,bi,bj)*e12p(i+1,j  ,bi,bj)*pi12(i+1,j  ,bi,bj)
     &   + rAz(i  ,j+1,bi,bj)*e12p(i  ,j+1,bi,bj)*pi12(i  ,j+1,bi,bj)
     &   + rAz(i+1,j+1,bi,bj)*e12p(i+1,j+1,bi,bj)*pi12(i+1,j+1,bi,bj) )
        denom = 1. / MAX(1.,SQRT(2.*(
     &       pi11(i,j,bi,bj)**2 + pi22(i,j,bi,bj)**2  + 2.*pi12C**2
     &       ))) / deltaRegSq
        c1(i,j) = ( dtau11*tau11 + dtau22*tau22
     &       + 2.*dtt12C ) * denom
CML     &       + 2.*dtau12C*tau12C ) * denom
        c2(i,j) = ( dtau11*pi11(i,j,bi,bj) + dtau22*pi22(i,j,bi,bj)
     &       + 2.*dtp12C ) * denom
CML     &       + 2.*dtau12C* pi12C ) * denom
        dpi11(i,j,bi,bj)= -pi11(i,j,bi,bj) + (tau11+dtau11)*recip_delta
     &       -  ( c1(i,j) * pi11(i,j,bi,bj) + c2(i,j)*tau11 )
        dpi22(i,j,bi,bj)= -pi22(i,j,bi,bj) + (tau22+dtau22)*recip_delta
     &       -  ( c1(i,j) * pi22(i,j,bi,bj) + c2(i,j)*tau22 )
        IF ( i.eq.20.and.j.eq.20 ) THEN !.and. myIter .le. 0) THEN
         print *, 'ml-myi', myIter
CML         print *, 'ml-de ', e11p(i,j,bi,bj), e22p(i,j,bi,bj)
CML     &        , e11(i,j,bi,bj), e22(i,j,bi,bj), c1(i,j), c2(i,j)
         print *, 'ml-dpi', dpi11(i,j,bi,bj), dpi22(i,j,bi,bj)
     &        , e11p(i,j,bi,bj), e22p(i,j,bi,bj), e12p(i,j,bi,bj)
CML     &        , denom, deltaReg
CML     &        , c1(i,j)*pi11(i,j,bi,bj), c2(i,j)*tau11
         print *, 'ml-pi ', pi11(i,j,bi,bj), pi22(i,j,bi,bj), pi12C
     &        , pi12(i,j,bi,bj)
         print *, 'ml-tau', tau11*recip_delta, tau22*recip_delta,
     &        tau12C*recip_delta, recip_ecc*e12(i,j,bi,bj)*recip_delta
         print *, 'ml-ttt', c1(i,j)*pi11(i,j,bi,bj), c2(i,j)*tau11,
     &        dtau11*recip_delta
         print *, 'ml-tt2', dtau11, dtau22, dtau12C, recip_delta
CML     & , dtau11*recip_delta, dtau22*recip_delta, dtau12C*recip_delta
CML         print *, 'ml-dp0', dpi11(i,j,bi,bj), pi11(i,j,bi,bj),
CML     &        tau11*recip_delta,dtau11*recip_delta
CML     &        , c1(i,j)*pi11(i,j,bi,bj), c2(i,j)*tau11
        ENDIF
       ENDDO
      ENDDO
      DO j=1,sNy+1
       DO i=1,sNx+1
C     some averaging to Z-points
        sumNorm  = HEFFM(i,j,  bi,bj)+HEFFM(i-1,j,  bi,bj)
     &           + HEFFM(i,j-1,bi,bj)+HEFFM(i-1,j-1,bi,bj)
        IF ( sumNorm.GT.0. _d 0 ) sumNorm = 1. _d 0 / sumNorm
#ifdef SEAICE_DELTA_SMOOTHREG
        deltaRegsq = sumNorm *
     &       ( deltaC(i,j  ,bi,bj)**2 + deltaC(i-1,j  ,bi,bj)**2
     &       + deltaC(i,j-1,bi,bj)**2 + deltaC(i-1,j-1,bi,bj)**2 )
        deltaRegSq = deltaRegSq + dsq
        deltaReg   = SQRT(deltaRegSq)
C        deltaReg = deltaC(i,j,bi,bj) + SEAICE_deltaMin
C        deltaRegSq = deltaReg**2
#else
        deltaReg = sumNorm *
     &       ( deltaC(i,j  ,bi,bj)  + deltaC(i-1,j  ,bi,bj)
     &       + deltaC(i,j-1,bi,bj)  + deltaC(i-1,j-1,bi,bj) )
        deltaReg   = MAX(deltaReg, SEAICE_deltaMin)
        deltaRegSq = deltaReg**2
#endif
        recip_delta = 1./deltaReg
        c1Z = sumNorm * (c1(i,j) + c1(i-1,j) + c1(i,j-1) + c1(i-1,j-1))
        c2Z = sumNorm * (c2(i,j) + c2(i-1,j) + c2(i,j-1) + c2(i-1,j-1))
        tau12  = recip_ecc*e12(i,j,bi,bj)
        dtau12 = recip_ecc*e12p(i,j,bi,bj)
        dpi12(i,j,bi,bj)= -pi12(i,j,bi,bj) + (tau12+dtau12)*recip_delta
     &       - ( c1Z * pi12(i,j,bi,bj) + c2Z * tau12 )
        IF ( i.eq.20.and.j.eq.20 .and. myIter .le. 0) THEN
         print *, 'ml-pi ', pi12(i,j,bi,bj)
         print *, 'ml-tau', tau12*recip_delta
        ENDIF
       ENDDO
      ENDDO
#endif /* SEAICE_CGRID */
      RETURN
      END
