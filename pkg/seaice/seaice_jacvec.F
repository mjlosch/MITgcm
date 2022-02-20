#include "SEAICE_OPTIONS.h"

CBOP
C     !ROUTINE: SEAICE_JACVEC
C     !INTERFACE:
      SUBROUTINE SEAICE_JACVEC(
     I     uIce, vIce, uIceRes, vIceRes,
     U     duIce, dvIce,
     I     newtonIter, krylovIter, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_JACVEC
C     | o For Jacobian-free Newton-Krylov solver compute
C     |   Jacobian times vector by finite difference approximation
C     *==========================================================*
C     | written by Martin Losch, Oct 2012
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"

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
C     u/vIce    :: local copies of the current ice velocity
      _RL uIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     u/vIceRes :: initial residual of this Newton iterate
      _RL uIceRes(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceRes(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     du/vIce   :: correction of ice velocities
      _RL duIce  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dvIce  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

#ifdef SEAICE_ALLOW_JFNK
C     Local variables:
      _RL utp     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vtp     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     u/vIceResP :: residual computed with u/vtp
      _RL uIceResP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceResP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C     i,j,bi,bj :: loop indices
      INTEGER i,j,bi,bj
      _RL epsilon, reps
CEOP
C     Instructions for using TAF or TAMC to generate exact Jacobian times
C     vector operations (if SEAICE_ALLOW_MOM_ADVECTION is defined, the
C     file list also needs to include seaice_mom_advection.f,
C     mom_calc_hfacz.f, mom_calc_ke.f, mom_calc_relvort3.f,
C     mom_vi_u_coriolis.f, mom_vi_u_coriolis_c4.f, mom_vi_u_grad_ke.f,
C     mom_vi_v_coriolis.f, mom_vi_v_coriolis_c4.f, mom_vi_v_grad_ke.f
C     plus flow information for diagnostics_fill.f:
CCCCCCCADJ SUBROUTINE DIAGNOSTICS_FILL INPUT  = 1,2,3,4,5,6,7,8
CCCCCCCADJ SUBROUTINE DIAGNOSTICS_FILL OUTPUT =
C     )
C
C     1. make small_f
C     2. cat seaice_calc_residual.f seaice_oceandrag_coeffs.f \
C        seaice_bottomdrag_coeffs.f seaice_calc_stressdiv.f \
C        seaice_calc_strainrates.f seaice_calc_viscosities.f \
C        seaice_calc_rhs.f seaice_calc_lhs.f > taf_input.f
C     3. staf -v1 -forward -toplevel seaice_calc_residual \
C             -input uIce,vice -output uIceRes,vIceRes taf_input.f
C     4. insert content of taf_input_ftl.f at the end of this file
C     5. add the following code and comment out the finite difference code
C
C     Instruction for using TAF 2.4 and higher (or staf with default -v2
C     starting with version 2.0):
C
C     1. make small_f
C     2. files="seaice_calc_residual.f seaice_oceandrag_coeffs.f \
C               seaice_bottomdrag_coeffs.f seaice_calc_stressdiv.f \
C               seaice_calc_strainrates.f seaice_calc_viscosities.f \
C               seaice_calc_rhs.f seaice_calc_lhs.f"
C     3. staf -forward -toplevel seaice_calc_residual \
C             -input uIce, vIce -output uIceRes,vIceRes $files
C     4. copy files seaice_*_tl.f to the corresponding seaice_*.f files,
C        e.g. with this bash script:
C     for file in $files; do
C       nfile=`echo $file | awk -F. '{printf "%s_tl.f", $1}'`;
C       \cp -f $nfile $file
C     done
C     5. add the following code, change "call g_seaice_calc_residual"
C        to "call seaice_calc_residual_tl", and comment out the finite
C        difference code
CML      _RL g_duIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL g_dvIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL g_uIceRes(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL g_vIceRes(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML
CMLC     Initialise
CML      DO bj=myByLo(myThid),myByHi(myThid)
CML       DO bi=myBxLo(myThid),myBxHi(myThid)
CML        DO J=1-Oly,sNy+Oly
CML         DO I=1-Olx,sNx+Olx
CML          g_duIce(I,J,bi,bj) = duice(I,J,bi,bj)
CML          g_dvIce(I,J,bi,bj) = dvice(I,J,bi,bj)
CML          g_uIceRes(I,J,bi,bj) = 0. _d 0
CML          g_vIceRes(I,J,bi,bj) = 0. _d 0
CML          uIceResP(I,J,bi,bj) = 0. _d 0
CML          vIceResP(I,J,bi,bj) = 0. _d 0
CML         ENDDO
CML        ENDDO
CML       ENDDO
CML      ENDDO
CML
CML      CALL G_SEAICE_CALC_RESIDUAL( uIce, g_duice, vIce,
CML     $g_dvice, uiceresp, g_uiceres, viceresp, g_viceres, newtoniter,
CML     $kryloviter, mytime, myiter, mythid )
CMLCML      For staf -v2 replace the above with the below call
CMLCML      CALL SEAICE_CALC_RESIDUAL_TL( uIce, g_duice, vIce,
CMLCML     $g_dvice, uiceresp, g_uiceres, viceresp, g_viceres, newtoniter,
CMLCML     $kryloviter, mytime, myiter, mythid )
CML
CML      DO bj=myByLo(myThid),myByHi(myThid)
CML       DO bi=myBxLo(myThid),myBxHi(myThid)
CML        DO J=1-Oly,sNy+Oly
CML         DO I=1-Olx,sNx+Olx
CML          duice(I,J,bi,bj)=g_uiceres(I,J,bi,bj)
CML          dvice(I,J,bi,bj)=g_viceres(I,J,bi,bj)
CML         ENDDO
CML        ENDDO
CML       ENDDO
CML      ENDDO

      IF ( SEAICEuseVelStress ) THEN
C     Compute Jacobian times vector following Shih et al.
       CALL SEAICE_JACVEC_VELSTRESS(
     I      uIce, vIce,
     U      duIce, dvIce,
     I      newtonIter, krylovIter, myTime, myIter, myThid )
      ELSE
C     Compute Jacobian times vector with pure finite difference
C     Initialise
       epsilon = SEAICE_JFNKepsilon
       reps    = 1. _d 0/epsilon

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO J=1-Oly,sNy+Oly
          DO I=1-Olx,sNx+Olx
           utp(I,J,bi,bj) = uIce(I,J,bi,bj) + epsilon * duIce(I,J,bi,bj)
           vtp(I,J,bi,bj) = vIce(I,J,bi,bj) + epsilon * dvIce(I,J,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO

C     Compute new residual F(u)
       CALL SEAICE_CALC_RESIDUAL(
     I      utp, vtp,
     O      uIceResP, vIceResP,
     I      newtonIter, krylovIter, myTime, myIter, myThid )

C     approximate Jacobian times vector by one-sided finite differences
C     and store in du/vIce
       DO bj = myByLo(myThid),myByHi(myThid)
        DO bi = myBxLo(myThid),myBxHi(myThid)
         DO I = 1, sNx
          DO J = 1, sNy
           duIce(I,J,bi,bj) =
     &          (uIceResP(I,J,bi,bj)-uIceRes(I,J,bi,bj))*reps
           dvIce(I,J,bi,bj) =
     &          (vIceResP(I,J,bi,bj)-vIceRes(I,J,bi,bj))*reps
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF

      RETURN
      END
C
      SUBROUTINE SEAICE_JACVEC_VELSTRESS(
     I     uIceLoc, vIceLoc,
     U     duIce, dvIce,
     I     newtonIter, krylovIter, myTime, myIter, myThid )
C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"
#include "DYNVARS.h"

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
C     du/vIce   :: velocity increment
      _RL duIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dvIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     Local variables:
      _RL uIceLHS (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceLHS (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL utp     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vtp     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e11p    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22p    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12p    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DWATNp  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12Csq  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL e12Csqp (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig11   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig22   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig11p  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig22p  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12p  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL strsDivX(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL strsDivY(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL zeros   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL coeff   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)

C     i,j,bi,bj :: loop indices
      INTEGER i,j,bi,bj
      INTEGER k
      _RL epsilon, reps
      _RS     SINWAT
      _RL     COSWAT, recip_deltaT
C     backward difference extrapolation factor
      _RL bdfAlpha
      _RL du   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL dv   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL dup  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL dvp  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C     fractional area at velocity points
      _RL areaW(1:sNx,1:sNy)
      _RL areaS(1:sNx,1:sNy)
#ifdef SEAICE_ALLOW_BOTTOMDRAG
C     perturbed basal drag coefficient
      _RL cBotCp  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
#ifdef SEAICE_ALLOW_MOM_ADVECTION
C     tendency due to advection of momentum
      _RL gUmom(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gVmom(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif /*  SEAICE_ALLOW_MOM_ADVECTION */
C     some helper variables
      _RL denominator, pip, piNorm
      _RL coeffZ, sumNorm
      _RL deltaRegSq, recip_ecc2
      _RL dragp, drag0
CEOP

C     Initialise
      epsilon = SEAICE_JFNKepsilon
      reps    = 1. _d 0/epsilon
      recip_ecc2 = 1./SEAICE_eccen**2

      k=1
      if ( usingPcoords ) k = Nr
      recip_deltaT = 1. _d 0 / SEAICE_deltaTdyn
C--   introduce turning angles
      SINWAT=SIN(SEAICE_waterTurnAngle*deg2rad)
      COSWAT=COS(SEAICE_waterTurnAngle*deg2rad)
C     backward difference extrapolation factor
      bdfAlpha = 1. _d 0
      IF ( SEAICEuseBDF2 ) THEN
       IF ( myIter.EQ.nIter0 .AND. SEAICEmomStartBDF.EQ.0 ) THEN
        bdfAlpha = 1. _d 0
       ELSE
        bdfAlpha = 1.5 _d 0
       ENDIF
      ENDIF

C     initialise fractional areas at velocity points
C     compute lhs side of momentum equations
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        IF ( SEAICEscaleSurfStress ) THEN
         DO j=1,sNy
          DO i=1,sNx
           areaW(i,j) = 0.5 _d 0*(AREA(i,j,bi,bj)+AREA(i-1,j,bi,bj))
           areaS(i,j) = 0.5 _d 0*(AREA(i,j,bi,bj)+AREA(i,j-1,bi,bj))
          ENDDO
         ENDDO
        ELSE
         DO j=1,sNy
          DO i=1,sNx
           areaW(i,j) = 1. _d 0
           areaS(i,j) = 1. _d 0
          ENDDO
         ENDDO
        ENDIF

        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          uIceLHS(i,j,bi,bj) = 0. _d 0
          vIceLHS(i,j,bi,bj) = 0. _d 0
          utp(i,j,bi,bj) = uIce(i,j,bi,bj) + epsilon * duIce(i,j,bi,bj)
          vtp(i,j,bi,bj) = vIce(i,j,bi,bj) + epsilon * dvIce(i,j,bi,bj)
          zeros(i,j,bi,bj) = 0.
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C     u/vIceLoc have changed so that new drag coefficients and
C     viscosities are required
      CALL SEAICE_OCEANDRAG_COEFFS(
     I     uIce, vIce, HEFFM,
     O     DWATN,
     I     krylovIter, myTime, myIter, myThid )
#ifdef SEAICE_ALLOW_BOTTOMDRAG
      CALL SEAICE_BOTTOMDRAG_COEFFS(
     I     uIce, vIce, HEFFM,
#ifdef SEAICE_ITD
     I     HEFFITD, AREAITD, AREA,
#else
     I     HEFF, AREA,
#endif
     O     CbotC,
     I     krylovIter, myTime, myIter, myThid )
#endif /* SEAICE_ALLOW_BOTTOMDRAG */
      CALL SEAICE_CALC_STRAINRATES(
     I     uIce, vIce,
     O     e11, e22, e12,
     I     krylovIter, myTime, myIter, myThid )
      CALL SEAICE_CALC_VISCOSITIES(
     I     e11, e22, e12, zMin, zMax, HEFFM, press0, tensileStrFac,
     O     eta, etaZ, zeta, zetaZ, press, deltaC,
     I     krylovIter, myTime, myIter, myThid )

C     repeat for some fields with u/vIce + epsilon * du/vIce
      CALL SEAICE_OCEANDRAG_COEFFS(
     I     utp, vtp, HEFFM,
     O     DWATNp,
     I     krylovIter, myTime, myIter, myThid )
#ifdef SEAICE_ALLOW_BOTTOMDRAG
      CALL SEAICE_BOTTOMDRAG_COEFFS(
     I     utp, vtp, HEFFM,
#ifdef SEAICE_ITD
     I     HEFFITD, AREAITD, AREA,
#else
     I     HEFF, AREA,
#endif
     O     CbotCp,
     I     krylovIter, myTime, myIter, myThid )
#endif /* SEAICE_ALLOW_BOTTOMDRAG */
      CALL SEAICE_CALC_STRAINRATES(
     I     duIce, dvIce,
     O     e11p, e22p, e12p,
     I     krylovIter, myTime, myIter, myThid )

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
C     stress components and linear perturbations with fixed coefficients
C     zeta and eta
        CALL SEAICE_CALC_STRESS(
     I       e11, e22, e12, zeros, zeta, eta, etaZ,
     O       sig11, sig22, sig12,
     I       bi, bj, myTime, myIter, myThid )
        CALL SEAICE_CALC_STRESS(
     I       e11p, e22p, e12p, zeros, zeta, eta, etaZ,
     O       sig11p, sig22p, sig12p,
     I       bi, bj, myTime, myIter, myThid )
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
C     area weighted average of the squares of e12 is more accurate
C     (and energy conserving)
          e12Csq(i,j) = 0.25 _d 0 * recip_rA(i,j,bi,bj) *
     &         ( rAz(i  ,j  ,bi,bj)* e12(i  ,j  ,bi,bj)**2
     &         + rAz(i+1,j  ,bi,bj)* e12(i+1,j  ,bi,bj)**2
     &         + rAz(i  ,j+1,bi,bj)* e12(i  ,j+1,bi,bj)**2
     &         + rAz(i+1,j+1,bi,bj)* e12(i+1,j+1,bi,bj)**2 )
          e12Csqp(i,j) = 0.25 _d 0 * recip_rA(i,j,bi,bj) *
     &         ( rAz(i  ,j  ,bi,bj)*e12p(i  ,j  ,bi,bj)**2
     &         + rAz(i+1,j  ,bi,bj)*e12p(i+1,j  ,bi,bj)**2
     &         + rAz(i  ,j+1,bi,bj)*e12p(i  ,j+1,bi,bj)**2
     &         + rAz(i+1,j+1,bi,bj)*e12p(i+1,j+1,bi,bj)**2 )
         ENDDO
        ENDDO
        DO j=0,sNy+1
         DO i=0,sNx+1
#ifdef SEAICE_DELTA_SMOOTHREG
          deltaRegSq  = (deltaC(i,j,bi,bj)**2 + SEAICE_deltaMin**2)
CML          deltaRegSq  = (deltaC(i,j,bi,bj)+ SEAICE_deltaMin)**2
#else
          deltaRegSq = MAX(deltaC(i,j,bi,bj),SEAICE_deltaMin)**2
#endif
          pip         = e11(i,j,bi,bj)+e22(i,j,bi,bj)
          piNorm      = (1.-recip_ecc2)*pip**2
     &         + 4.*e12Csq(i,j)*recip_ecc2
          piNorm      = 2. _d 0 * SQRT(piNorm/deltaRegSq)
CML          print *, 'ml-piNorm =', piNorm
          denominator = MAX(1.D0,piNorm) * deltaRegSq
C     this is the eps:tilde(eps) part
          coeff(i,j)  = e11(i,j,bi,bj)*e11p(i,j,bi,bj)
     &         +        e22(i,j,bi,bj)*e22p(i,j,bi,bj)
     &         + 2. _d 0 * SQRT(e12Csq(i,j)*e12Csqp(i,j))
C     this is the tr(eps):tr(tilde(eps)) part
          coeff(i,j) = ( 2. _d 0*coeff(i,j)*recip_ecc2
     &         + (  e11(i,j,bi,bj)+ e22(i,j,bi,bj) )
     &         * ( e11p(i,j,bi,bj)+e22p(i,j,bi,bj) )
     &         )/denominator
         ENDDO
        ENDDO
        DO j=0,sNy
         DO i=0,sNx
          sig11p(i,j) = sig11p(i,j) - coeff(i,j)*sig11(i,j)
          sig22p(i,j) = sig22p(i,j) - coeff(i,j)*sig22(i,j)
         ENDDO
        ENDDO
        DO j=1,sNy+1
         DO i=1,sNx+1
          sumNorm  = HEFFM(i,j,  bi,bj)+HEFFM(i-1,j,  bi,bj)
     &             + HEFFM(i,j-1,bi,bj)+HEFFM(i-1,j-1,bi,bj)
          IF ( sumNorm.GT.0. _d 0 ) sumNorm = 1. _d 0 / sumNorm
          coeffz = sumNorm*( coeff(i,j)  +coeff(i-1,j)
     &                     + coeff(i,j-1)+coeff(i-1,j-1) )
          sig12p(i,j) = sig12p(i,j) - coeffz*sig12(i,j)
         ENDDO
        ENDDO
C     compute "stress divergence"
        DO j=1,sNy
         DO i=1,sNx
          strsDivX(i,j) =
     &         ( sig11p(i  ,j  ) * _dyF(i  ,j,bi,bj)
     &         - sig11p(i-1,j  ) * _dyF(i-1,j,bi,bj)
     &         + sig12p(i  ,j+1) * _dxV(i,j+1,bi,bj)
     &         - sig12p(i  ,j  ) * _dxV(i,j  ,bi,bj)
     &         ) * recip_rAw(i,j,bi,bj)
          strsDivY(i,j) =
     &         ( sig22p(i  ,j  ) * _dxF(i,j  ,bi,bj)
     &         - sig22p(i  ,j-1) * _dxF(i,j-1,bi,bj)
     &         + sig12p(i+1,j  ) * _dyU(i+1,j,bi,bj)
     &         - sig12p(i  ,j  ) * _dyU(i  ,j,bi,bj)
     &         ) * recip_rAs(i,j,bi,bj)
         ENDDO
        ENDDO
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          dup(i,j) = utp (i,j,bi,bj) - uVel(i,j,k,bi,bj)
          dvp(i,j) = vtp (i,j,bi,bj) - vVel(i,j,k,bi,bj)
          du (i,j) = uIce(i,j,bi,bj) - uVel(i,j,k,bi,bj)
          dv (i,j) = vIce(i,j,bi,bj) - vVel(i,j,k,bi,bj)
         ENDDO
        ENDDO
        DO j=1,sNy
         DO i=1,sNx
C     mass*(uIce)/deltaT - dsigma/dx
          uIceLHS(i,j,bi,bj) =
     &         bdfAlpha*seaiceMassU(i,j,bi,bj)*recip_deltaT
     &         *duIce(i,j,bi,bj) - strsDivX(i,j)
C     mass*(vIce)/deltaT - dsigma/dy
          vIceLHS(i,j,bi,bj) =
     &         bdfAlpha*seaiceMassV(i,j,bi,bj)*recip_deltaT
     &         *dvIce(i,j,bi,bj) - strsDivY(i,j)
C     coriols terms: - mass*f*vIce
          uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj) - 0.5 _d 0*(
     &         seaiceMassC(i  ,j,bi,bj) * _fCori(i  ,j,bi,bj)
     &       * 0.5 _d 0*( dvIce(i  ,j,bi,bj)+dvIce(i  ,j+1,bi,bj) )
     &       + seaiceMassC(i-1,j,bi,bj) * _fCori(i-1,j,bi,bj)
     &       * 0.5 _d 0*( dvIce(i-1,j,bi,bj)+dvIce(i-1,j+1,bi,bj) )
     &           )
C                    + mass*f*uIce
          vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj) + 0.5 _d 0*(
     &         seaiceMassC(i,j  ,bi,bj) * _fCori(i,j  ,bi,bj)
     &       * 0.5 _d 0*( duIce(i,j  ,bi,bj)+duIce(i+1,  j,bi,bj) )
     &       + seaiceMassC(i,j-1,bi,bj) * _fCori(i,j-1,bi,bj)
     &       * 0.5 _d 0*( duIce(i,j-1,bi,bj)+duIce(i+1,j-1,bi,bj) )
     &           )
C     ocean-ice and bottom drag terms: + (Cdrag*cosWat+Cb)*uIce - vIce*sinWat)
          dragp = (
#ifdef SEAICE_ALLOW_BOTTOMDRAG
     &         0.5 _d 0 * ( CbotCp(i,j,bi,bj)+CbotCp(i-1,j,bi,bj) )
     &         * utp(i,j,bi,bj) +
#endif
     &         0.5 _d 0 * ( DWATNp(i,j,bi,bj)+DWATNp(i-1,j,bi,bj) )
     &         * COSWAT * dup(i,j)
     &         - SIGN(SINWAT, _fCori(i,j,bi,bj))* 0.5 _d 0 *
     &         ( DWATNp(i  ,j,bi,bj) *0.5 _d 0*(dvp(i  ,j)+dvp(i  ,j+1))
     &         + DWATNp(i-1,j,bi,bj) *0.5 _d 0*(dvp(i-1,j)+dvp(i-1,j+1))
     &         ) ) * areaW(i,j)
          drag0 = (
#ifdef SEAICE_ALLOW_BOTTOMDRAG
     &         0.5 _d 0 * ( CbotC(i,j,bi,bj)+CbotC(i-1,j,bi,bj) )
     &         * uIce(i,j,bi,bj) +
#endif
     &         0.5 _d 0 * ( DWATN(i,j,bi,bj)+DWATN(i-1,j,bi,bj) )
     &         * COSWAT * du(i,j)
     &         - SIGN(SINWAT, _fCori(i,j,bi,bj))* 0.5 _d 0 *
     &         ( DWATN(i  ,j,bi,bj) * 0.5 _d 0 * (dv(i  ,j)+dv(i  ,j+1))
     &         + DWATN(i-1,j,bi,bj) * 0.5 _d 0 * (dv(i-1,j)+dv(i-1,j+1))
     &         ) ) * areaW(i,j)
          uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj) + (dragp-drag0)*reps
C                                      + (Cdrag*cosWat+Cb)*uIce + uIce*sinWat)
          dragp = (
#ifdef SEAICE_ALLOW_BOTTOMDRAG
     &         0.5 _d 0 * ( CbotCp(i,j,bi,bj)+CbotCp(i,j-1,bi,bj) )
     &         * vtp(i,j,bi,bj) +
#endif
     &         0.5 _d 0 * ( DWATNp(i,j,bi,bj)+DWATNp(i,j-1,bi,bj) )
     &         * COSWAT * dvp(i,j)
     &         + SIGN(SINWAT, _fCori(i,j,bi,bj)) * 0.5 _d 0 *
     &         ( DWATNp(i,j  ,bi,bj) *0.5 _d 0*(dup(i,j  )+dup(i+1,j  ))
     &         + DWATNp(i,j-1,bi,bj) *0.5 _d 0*(dup(i,j-1)+dup(i+1,j-1))
     &         ) ) * areaS(i,j)
          drag0 = (
#ifdef SEAICE_ALLOW_BOTTOMDRAG
     &         0.5 _d 0 * ( CbotC(i,j,bi,bj)+CbotC(i,j-1,bi,bj) )
     &         * vIce(i,j,bi,bj) +
#endif
     &         0.5 _d 0 * ( DWATN(i,j,bi,bj)+DWATN(i,j-1,bi,bj) )
     &         * COSWAT * dv(i,j)
     &         + SIGN(SINWAT, _fCori(i,j,bi,bj)) * 0.5 _d 0 *
     &         ( DWATN(i,j  ,bi,bj) * 0.5 _d 0 * (du(i,j  )+du(i+1,j  ))
     &         + DWATN(i,j-1,bi,bj) * 0.5 _d 0 * (du(i,j-1)+du(i+1,j-1))
     &         ) ) * areaS(i,j)
          vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj) + (dragp-drag0)*reps
C     apply masks for interior (important when we have open boundaries)
          uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj)*maskinW(i,j,bi,bj)
          vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj)*maskinS(i,j,bi,bj)
         ENDDO
        ENDDO
#ifdef SEAICE_ALLOW_MOM_ADVECTION
        IF ( SEAICEmomAdvection ) THEN
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           gUmom(i,j) = 0. _d 0
           gVmom(i,j) = 0. _d 0
          ENDDO
         ENDDO
         CALL SEAICE_MOM_ADVECTION(
     I        bi,bj,1,sNx,1,sNy,
     I        utp, vtp,
     O        gUmom, gVmom,
     I        myTime, myIter, myThid )
C     Beware of sign! gU/Vmom is computed for the rhs of the equation;
C     therefore, we need to substract gU/Vmom from the left hand side
         DO j=1,sNy
          DO i=1,sNx
           uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj) - gUmom(i,j)*reps
           vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj) - gVmom(i,j)*reps
          ENDDO
         ENDDO
         CALL SEAICE_MOM_ADVECTION(
     I        bi,bj,1,sNx,1,sNy,
     I        uIceLoc, vIceLoc,
     O        gUmom, gVmom,
     I        myTime, myIter, myThid )
         DO j=1,sNy
          DO i=1,sNx
           uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj) + gUmom(i,j)*reps
           vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj) + gVmom(i,j)*reps
          ENDDO
         ENDDO
        ENDIF
#endif /* SEAICE_ALLOW_MOM_ADVECTION */
        DO j=1,sNy
         DO i=1,sNx
           duIce(i,j,bi,bj) = uIceLHS(i,j,bi,bj)
           dvIce(i,j,bi,bj) = vIceLHS(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#endif /* SEAICE_ALLOW_JFNK */

      RETURN
      END
