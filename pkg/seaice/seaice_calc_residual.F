#include "SEAICE_OPTIONS.h"

CBOP
C     !ROUTINE: SEAICE_CALC_RESIDUAL
C     !INTERFACE:
      SUBROUTINE SEAICE_CALC_RESIDUAL(
     I     uIceLoc, vIceLoc,
     O     uIceRes, vIceRes,
     I     newtonIter, krylovIter, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_CALC_RESIDUAL
C     | o For Jacobian-free Newton-Krylov solver compute
C     |   the residual of the momentum equations
C     *==========================================================*
C     | written by Martin Losch, Oct 2012
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
C     newtonIter :: current iterate of Newton iteration
C     krylovIter :: current iterate of Krylov iteration
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER newtonIter
      INTEGER krylovIter
C     u/vIceLoc :: local copies of the current ice velocity
      _RL uIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     u/vIceRes :: residual of sea-ice momentum equations
      _RL uIceRes(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceRes(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

#ifdef SEAICE_ALLOW_JFNK
C     u/vIceLHS :: left hand side of momentum equations
      _RL uIceLHS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceLHS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     u/vIceRHS :: righ hand side of momentum equations
      _RL uIceRHS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceRHS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C     i,j,bi,bj :: loop indices
      INTEGER i,j,bi,bj
CEOP

C     Initialise
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO J=1-OLy,sNy+OLy
         DO I=1-OLx,sNx+OLx
          uIceLHS(I,J,bi,bj) = 0. _d 0
          vIceLHS(I,J,bi,bj) = 0. _d 0
          uIceRHS(I,J,bi,bj) = 0. _d 0
          vIceRHS(I,J,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C     u/vIceLoc have changed so that new drag coefficients and
C     viscosities are required
      CALL SEAICE_OCEANDRAG_COEFFS(
     I     uIceLoc, vIceLoc, HEFFM,
     O     DWATN,
     I     krylovIter, myTime, myIter, myThid )
#ifdef SEAICE_ALLOW_BOTTOMDRAG
      CALL SEAICE_BOTTOMDRAG_COEFFS(
     I     uIceLoc, vIceLoc, HEFFM,
#ifdef SEAICE_ITD
     I     HEFFITD, AREAITD, AREA,
#else
     I     HEFF, AREA,
#endif
     O     CbotC,
     I     krylovIter, myTime, myIter, myThid )
#endif /* SEAICE_ALLOW_BOTTOMDRAG */
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO J=1-Oly,sNy+Oly
         DO I=1-Olx,sNx+Olx
          e11p(I,J,bi,bj) = 0. _d 0
          e22p(I,J,bi,bj) = 0. _d 0
          e12p(I,J,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      IF ( SEAICEuseVelStress .AND. krylovIter .GE. 0 ) THEN
C     compute strain rates for unperturbed velocities
       CALL SEAICE_CALC_STRAINRATES(
     I      uIce, vIce,
     O      e11, e22, e12,
     I      krylovIter, myTime, myIter, myThid )
C     and for the perturbations velocities (stored in u/vIceRes)
       CALL SEAICE_CALC_STRAINRATES(
     I      uIceRes, vIceRes,
     O      e11p, e22p, e12p,
     I      krylovIter, myTime, myIter, myThid )
      ELSE
       CALL SEAICE_CALC_STRAINRATES(
     I      uIceLoc, vIceLoc,
     O      e11, e22, e12,
     I      krylovIter, myTime, myIter, myThid )
      ENDIF
      CALL SEAICE_CALC_VISCOSITIES(
     I     e11, e22, e12, zMin, zMax, HEFFM, press0, tensileStrFac,
     O     eta, etaZ, zeta, zetaZ, press, deltaC,
     I     krylovIter, myTime, myIter, myThid )

C     The scheme is backward Euler in time, i.e. the rhs-vector contains
C     only terms that are independent of u/vIce, except for the time
C     derivative part mass*(u/vIce-u/vIceNm1)/deltaT

C     compute new right hand side (depends to DWATN=Cdrag)
C     sea-surface tilt and wind stress: FORCEX0, FORCEY0
C     + mass*(u/vIceNm1)/deltaT
C     + Cdrag*(uVel*cosWat - vVel*sinWat)/(vVel*cosWat + uVel*sinWat)
      CALL SEAICE_CALC_RHS(
     O      uIceRHS, vIceRHS,
     I      newtonIter, krylovIter, myTime, myIter, myThid )

C     Left-hand side contributions:
C     + mass*(u/vIce)/deltaT
C     + Cdrag*(uIce*cosWat - vIce*sinWat)/(vIce*cosWat + uIce*sinWat)
C     + CdragBot*uIce/vIce
C     - mass*f*vIce/+mass*f*uIce
C     - dsigma/dx / -dsigma/dy, eta and zeta are only computed once per
C     Newton iterate
      IF ( SEAICEuseVelStress .AND. krylovIter .GE. 0 ) THEN
C     this is a hack for now to avoid passing more arguments to
C     seaice_calc_residual, because in this case u/vIceRes = du/vIce
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO J=1-Oly,sNy+Oly
          DO I=1-Olx,sNx+Olx
           uIceLHS(I,J,bi,bj) = uIceRes(I,J,bi,bj)
           vIceLHS(I,J,bi,bj) = vIceRes(I,J,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
CML      ELSE
CML       DO bj=myByLo(myThid),myByHi(myThid)
CML        DO bi=myBxLo(myThid),myBxHi(myThid)
CML         DO J=1-Oly,sNy+Oly
CML          DO I=1-Olx,sNx+Olx
CML           uIceLHS(I,J,bi,bj) = 0. _d 0
CML           vIceLHS(I,J,bi,bj) = 0. _d 0
CML          ENDDO
CML         ENDDO
CML        ENDDO
CML       ENDDO
      ENDIF
      CALL SEAICE_CALC_LHS(
     I     uIceLoc, vIceLoc,
     O     uIceLHS, vIceLHS,
     I     newtonIter, krylovIter, myTime, myIter, myThid )


C     Calculate the residual
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO J=1,sNy
          DO I=1,sNx
           uIceRes(I,J,bi,bj) = uIceLHS(I,J,bi,bj) - uIceRHS(I,J,bi,bj)
           vIceRes(I,J,bi,bj) = vIceLHS(I,J,bi,bj) - vIceRHS(I,J,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO

#endif /* SEAICE_ALLOW_JFNK */

      RETURN
      END
