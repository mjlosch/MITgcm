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
CMLC     Compute Jacobian times vector following Shih et al.
CML       CALL SEAICE_JACVEC_VELSTRESS(
CML     I      uIce, vIce,
CML     U      duIce, dvIce,
CML     I      newtonIter, krylovIter, myTime, myIter, myThid )
C     this is a hack for now to avoid passing more arguments to
C     seaice_calc_residual
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO J=1-Oly,sNy+Oly
          DO I=1-Olx,sNx+Olx
           uIceResP(I,J,bi,bj) = duIce(I,J,bi,bj)
           vIceResP(I,J,bi,bj) = dvIce(I,J,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF
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
     I     utp, vtp,
     O     uIceResP, vIceResP,
     I     newtonIter, krylovIter, myTime, myIter, myThid )

C     approximate Jacobian times vector by one-sided finite differences
C     and store in du/vIce
      DO bj = myByLo(myThid),myByHi(myThid)
       DO bi = myBxLo(myThid),myBxHi(myThid)
        DO I = 1, sNx
         DO J = 1, sNy
          duIce(I,J,bi,bj) =
     &         (uIceResP(I,J,bi,bj)-uIceRes(I,J,bi,bj))*reps
          dvIce(I,J,bi,bj) =
     &         (vIceResP(I,J,bi,bj)-vIceRes(I,J,bi,bj))*reps
         ENDDO
        ENDDO
       ENDDO
      ENDDO

CML      RETURN
CML      END
CMLC
CML      SUBROUTINE SEAICE_JACVEC_VELSTRESS(
CML     I     uIceLoc, vIceLoc,
CML     U     duIce, dvIce,
CML     I     newtonIter, krylovIter, myTime, myIter, myThid )
CMLC     !USES:
CML      IMPLICIT NONE
CML
CMLC     === Global variables ===
CML#include "SIZE.h"
CML#include "EEPARAMS.h"
CML#include "PARAMS.h"
CML#include "GRID.h"
CML#include "SEAICE_SIZE.h"
CML#include "SEAICE_PARAMS.h"
CML#include "SEAICE.h"
CML#include "DYNVARS.h"
CML
CMLC     !INPUT/OUTPUT PARAMETERS:
CMLC     === Routine arguments ===
CMLC     myTime :: Simulation time
CMLC     myIter :: Simulation timestep number
CMLC     myThid :: my Thread Id. number
CMLC     newtonIter :: current iterate of Newton iteration
CMLC     krylovIter :: current iterate of Krylov iteration
CML      _RL     myTime
CML      INTEGER myIter
CML      INTEGER myThid
CML      INTEGER newtonIter
CML      INTEGER krylovIter
CMLC     u/vIceLoc :: local copies of the current ice velocity
CML      _RL uIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL vIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CMLC     du/vIce   :: velocity increment
CML      _RL duIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL dvIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CMLC     Local variables:
CML      _RL uIceLHS (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL vIceLHS (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL utp     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL vtp     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL e11p    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL e22p    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL e12p    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL DWATNp  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL e12Csq  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL e12Csqp (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL sig11   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL sig22   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL sig12   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL sig11p  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL sig22p  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL sig12p  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL strsDivX(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL strsDivY(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL zeros   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML      _RL coeff   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML
CMLC     i,j,bi,bj :: loop indices
CML      INTEGER i,j,bi,bj
CML      INTEGER k
CML      _RL epsilon, reps
CML      _RS     SINWAT
CML      _RL     COSWAT, recip_deltaT
CMLC     backward difference extrapolation factor
CML      _RL bdfAlpha
CML      _RL du   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL dv   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL dup  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL dvp  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CMLC     fractional area at velocity points
CML      _RL areaW(1:sNx,1:sNy)
CML      _RL areaS(1:sNx,1:sNy)
CML#ifdef SEAICE_ALLOW_BOTTOMDRAG
CMLC     perturbed basal drag coefficient
CML      _RL cBotCp  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CML#endif
CML#ifdef SEAICE_ALLOW_MOM_ADVECTION
CMLC     tendency due to advection of momentum
CML      _RL gUmom(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML      _RL gVmom(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CML#endif /*  SEAICE_ALLOW_MOM_ADVECTION */
CMLC     some helper variables
CML      _RL denominator, pip, piNorm
CML      _RL coeffZ, sumNorm
CML      _RL deltaRegSq, recip_ecc2
CML      _RL dragp, drag0
CMLCEOP
CML
CMLC     Initialise
CML      epsilon = SEAICE_JFNKepsilon
CML      reps    = 1. _d 0/epsilon
CML      recip_ecc2 = 1./SEAICE_eccen**2
CML
CML      k=1
CML      if ( usingPcoords ) k = Nr
CML      recip_deltaT = 1. _d 0 / SEAICE_deltaTdyn
CMLC--   introduce turning angles
CML      SINWAT=SIN(SEAICE_waterTurnAngle*deg2rad)
CML      COSWAT=COS(SEAICE_waterTurnAngle*deg2rad)
CMLC     backward difference extrapolation factor
CML      bdfAlpha = 1. _d 0
CML      IF ( SEAICEuseBDF2 ) THEN
CML       IF ( myIter.EQ.nIter0 .AND. SEAICEmomStartBDF.EQ.0 ) THEN
CML        bdfAlpha = 1. _d 0
CML       ELSE
CML        bdfAlpha = 1.5 _d 0
CML       ENDIF
CML      ENDIF
CML
CMLC     initialise fractional areas at velocity points
CMLC     compute lhs side of momentum equations
CML      DO bj=myByLo(myThid),myByHi(myThid)
CML       DO bi=myBxLo(myThid),myBxHi(myThid)
CML        IF ( SEAICEscaleSurfStress ) THEN
CML         DO j=1,sNy
CML          DO i=1,sNx
CML           areaW(i,j) = 0.5 _d 0*(AREA(i,j,bi,bj)+AREA(i-1,j,bi,bj))
CML           areaS(i,j) = 0.5 _d 0*(AREA(i,j,bi,bj)+AREA(i,j-1,bi,bj))
CML          ENDDO
CML         ENDDO
CML        ELSE
CML         DO j=1,sNy
CML          DO i=1,sNx
CML           areaW(i,j) = 1. _d 0
CML           areaS(i,j) = 1. _d 0
CML          ENDDO
CML         ENDDO
CML        ENDIF
CML
CML        DO j=1-OLy,sNy+OLy
CML         DO i=1-OLx,sNx+OLx
CML          uIceLHS(i,j,bi,bj) = 0. _d 0
CML          vIceLHS(i,j,bi,bj) = 0. _d 0
CML          utp(i,j,bi,bj) = uIceLoc(i,j,bi,bj) + epsilon*duIce(i,j,bi,bj)
CML          vtp(i,j,bi,bj) = vIceLoc(i,j,bi,bj) + epsilon*dvIce(i,j,bi,bj)
CML          zeros(i,j,bi,bj) = 0.
CML         ENDDO
CML        ENDDO
CML       ENDDO
CML      ENDDO
CML
CMLC     u/vIceLoc have changed so that new drag coefficients and
CMLC     viscosities are required
CML      CALL SEAICE_OCEANDRAG_COEFFS(
CML     I     uIceLoc, vIceLoc, HEFFM,
CML     O     DWATN,
CML     I     krylovIter, myTime, myIter, myThid )
CML#ifdef SEAICE_ALLOW_BOTTOMDRAG
CML      CALL SEAICE_BOTTOMDRAG_COEFFS(
CML     I     uIceLoc, vIceLoc, HEFFM,
CML#ifdef SEAICE_ITD
CML     I     HEFFITD, AREAITD, AREA,
CML#else
CML     I     HEFF, AREA,
CML#endif
CML     O     CbotC,
CML     I     krylovIter, myTime, myIter, myThid )
CML#endif /* SEAICE_ALLOW_BOTTOMDRAG */
CML      CALL SEAICE_CALC_STRAINRATES(
CML     I     uIceLoc, vIceLoc,
CML     O     e11, e22, e12,
CML     I     krylovIter, myTime, myIter, myThid )
CML      CALL SEAICE_CALC_VISCOSITIES(
CML     I     e11, e22, e12, zMin, zMax, HEFFM, press0, tensileStrFac,
CML     O     eta, etaZ, zeta, zetaZ, press, deltaC,
CML     I     krylovIter, myTime, myIter, myThid )
CML
CMLC     repeat for some fields with u/vIce + epsilon * du/vIce
CML      CALL SEAICE_OCEANDRAG_COEFFS(
CML     I     utp, vtp, HEFFM,
CML     O     DWATNp,
CML     I     krylovIter, myTime, myIter, myThid )
CML#ifdef SEAICE_ALLOW_BOTTOMDRAG
CML      CALL SEAICE_BOTTOMDRAG_COEFFS(
CML     I     utp, vtp, HEFFM,
CML#ifdef SEAICE_ITD
CML     I     HEFFITD, AREAITD, AREA,
CML#else
CML     I     HEFF, AREA,
CML#endif
CML     O     CbotCp,
CML     I     krylovIter, myTime, myIter, myThid )
CML#endif /* SEAICE_ALLOW_BOTTOMDRAG */
CML      CALL SEAICE_CALC_STRAINRATES(
CML     I     duIce, dvIce,
CML     O     e11p, e22p, e12p,
CML     I     krylovIter, myTime, myIter, myThid )
CML
CML      DO bj=myByLo(myThid),myByHi(myThid)
CML       DO bi=myBxLo(myThid),myBxHi(myThid)
CMLC     stress components and linear perturbations with fixed coefficients
CMLC     zeta and eta
CML        CALL SEAICE_CALC_STRESS(
CML     I       e11, e22, e12, zeros, zeta, eta, etaZ,
CML     O       sig11, sig22, sig12,
CML     I       bi, bj, myTime, myIter, myThid )
CML        CALL SEAICE_CALC_STRESS(
CML     I       e11p, e22p, e12p, zeros, zeta, eta, etaZ,
CML     O       sig11p, sig22p, sig12p,
CML     I       bi, bj, myTime, myIter, myThid )
CML        DO j=1-OLy+1,sNy+OLy-1
CML         DO i=1-OLx+1,sNx+OLx-1
CMLC     area weighted average of the squares of e12 is more accurate
CMLC     (and energy conserving)
CML          e12Csq(i,j) = 0.25 _d 0 * recip_rA(i,j,bi,bj) *
CML     &         ( rAz(i  ,j  ,bi,bj)* e12(i  ,j  ,bi,bj)**2
CML     &         + rAz(i+1,j  ,bi,bj)* e12(i+1,j  ,bi,bj)**2
CML     &         + rAz(i  ,j+1,bi,bj)* e12(i  ,j+1,bi,bj)**2
CML     &         + rAz(i+1,j+1,bi,bj)* e12(i+1,j+1,bi,bj)**2 )
CML          e12Csqp(i,j) = 0.25 _d 0 * recip_rA(i,j,bi,bj) *
CML     &         ( rAz(i  ,j  ,bi,bj)*e12p(i  ,j  ,bi,bj)**2
CML     &         + rAz(i+1,j  ,bi,bj)*e12p(i+1,j  ,bi,bj)**2
CML     &         + rAz(i  ,j+1,bi,bj)*e12p(i  ,j+1,bi,bj)**2
CML     &         + rAz(i+1,j+1,bi,bj)*e12p(i+1,j+1,bi,bj)**2 )
CML         ENDDO
CML        ENDDO
CML        DO j=0,sNy+1
CML         DO i=0,sNx+1
CML#ifdef SEAICE_DELTA_SMOOTHREG
CML          deltaRegSq  = (deltaC(i,j,bi,bj)**2 + SEAICE_deltaMin**2)
CMLCML          deltaRegSq  = (deltaC(i,j,bi,bj)+ SEAICE_deltaMin)**2
CML#else
CML          deltaRegSq = MAX(deltaC(i,j,bi,bj),SEAICE_deltaMin)**2
CML#endif
CML          pip         = e11(i,j,bi,bj)+e22(i,j,bi,bj)
CML          piNorm      = (1.-recip_ecc2)*pip**2
CML     &         + 4.*e12Csq(i,j)*recip_ecc2
CML          piNorm      = 2. _d 0 * SQRT(piNorm/deltaRegSq)
CMLCML          print *, 'ml-piNorm =', piNorm
CMLCML          denominator = MAX(1.D0,piNorm) * deltaRegSq
CML          denominator = deltaRegSq
CMLC     this is the eps:tilde(eps) part
CML          coeff(i,j)  = e11(i,j,bi,bj)*e11p(i,j,bi,bj)
CML     &         +        e22(i,j,bi,bj)*e22p(i,j,bi,bj)
CML     &         + 2. _d 0 * SQRT(e12Csq(i,j)*e12Csqp(i,j))
CMLC     this is the tr(eps):tr(tilde(eps)) part
CML          coeff(i,j) = ( 2. _d 0*coeff(i,j)*recip_ecc2
CML     &         + ( 1. _d 0 - recip_ecc2 )
CML     &         * (  e11(i,j,bi,bj)+ e22(i,j,bi,bj) )
CML     &         * ( e11p(i,j,bi,bj)+e22p(i,j,bi,bj) )
CML     &         )/denominator
CML         ENDDO
CML        ENDDO
CML        DO j=0,sNy
CML         DO i=0,sNx
CML          sig11p(i,j) = sig11p(i,j) - coeff(i,j)*sig11(i,j)
CML          sig22p(i,j) = sig22p(i,j) - coeff(i,j)*sig22(i,j)
CML         ENDDO
CML        ENDDO
CML        DO j=1,sNy+1
CML         DO i=1,sNx+1
CML          sumNorm  = HEFFM(i,j,  bi,bj)+HEFFM(i-1,j,  bi,bj)
CML     &             + HEFFM(i,j-1,bi,bj)+HEFFM(i-1,j-1,bi,bj)
CML          IF ( sumNorm.GT.0. _d 0 ) sumNorm = 1. _d 0 / sumNorm
CML          coeffz = sumNorm*( coeff(i,j)  +coeff(i-1,j)
CML     &                     + coeff(i,j-1)+coeff(i-1,j-1) )
CML          sig12p(i,j) = sig12p(i,j) - coeffz*sig12(i,j)
CML         ENDDO
CML        ENDDO
CMLC     compute "stress divergence"
CML        DO j=1,sNy
CML         DO i=1,sNx
CML          strsDivX(i,j) =
CML     &         ( sig11p(i  ,j  ) * _dyF(i  ,j,bi,bj)
CML     &         - sig11p(i-1,j  ) * _dyF(i-1,j,bi,bj)
CML     &         + sig12p(i  ,j+1) * _dxV(i,j+1,bi,bj)
CML     &         - sig12p(i  ,j  ) * _dxV(i,j  ,bi,bj)
CML     &         ) * recip_rAw(i,j,bi,bj)
CML          strsDivY(i,j) =
CML     &         ( sig22p(i  ,j  ) * _dxF(i,j  ,bi,bj)
CML     &         - sig22p(i  ,j-1) * _dxF(i,j-1,bi,bj)
CML     &         + sig12p(i+1,j  ) * _dyU(i+1,j,bi,bj)
CML     &         - sig12p(i  ,j  ) * _dyU(i  ,j,bi,bj)
CML     &         ) * recip_rAs(i,j,bi,bj)
CML         ENDDO
CML        ENDDO
CML        DO j=1-OLy,sNy+OLy
CML         DO i=1-OLx,sNx+OLx
CML          dup(i,j) = utp (i,j,bi,bj)    - uVel(i,j,k,bi,bj)
CML          dvp(i,j) = vtp (i,j,bi,bj)    - vVel(i,j,k,bi,bj)
CML          du (i,j) = uIceLoc(i,j,bi,bj) - uVel(i,j,k,bi,bj)
CML          dv (i,j) = vIceLoc(i,j,bi,bj) - vVel(i,j,k,bi,bj)
CML         ENDDO
CML        ENDDO
CML        DO j=1,sNy
CML         DO i=1,sNx
CMLC     mass*(uIce)/deltaT - dsigma/dx
CML          uIceLHS(i,j,bi,bj) =
CML     &         bdfAlpha*seaiceMassU(i,j,bi,bj)*recip_deltaT
CML     &         *duIce(i,j,bi,bj) - strsDivX(i,j)
CMLC     mass*(vIce)/deltaT - dsigma/dy
CML          vIceLHS(i,j,bi,bj) =
CML     &         bdfAlpha*seaiceMassV(i,j,bi,bj)*recip_deltaT
CML     &         *dvIce(i,j,bi,bj) - strsDivY(i,j)
CMLC     coriols terms: - mass*f*vIce
CML          uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj) - 0.5 _d 0*(
CML     &         seaiceMassC(i  ,j,bi,bj) * _fCori(i  ,j,bi,bj)
CML     &       * 0.5 _d 0*( dvIce(i  ,j,bi,bj)+dvIce(i  ,j+1,bi,bj) )
CML     &       + seaiceMassC(i-1,j,bi,bj) * _fCori(i-1,j,bi,bj)
CML     &       * 0.5 _d 0*( dvIce(i-1,j,bi,bj)+dvIce(i-1,j+1,bi,bj) )
CML     &           )
CMLC                    + mass*f*uIce
CML          vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj) + 0.5 _d 0*(
CML     &         seaiceMassC(i,j  ,bi,bj) * _fCori(i,j  ,bi,bj)
CML     &       * 0.5 _d 0*( duIce(i,j  ,bi,bj)+duIce(i+1,  j,bi,bj) )
CML     &       + seaiceMassC(i,j-1,bi,bj) * _fCori(i,j-1,bi,bj)
CML     &       * 0.5 _d 0*( duIce(i,j-1,bi,bj)+duIce(i+1,j-1,bi,bj) )
CML     &           )
CMLC     ocean-ice and bottom drag terms: + (Cdrag*cosWat+Cb)*uIce - vIce*sinWat)
CML          dragp = (
CML#ifdef SEAICE_ALLOW_BOTTOMDRAG
CML     &         0.5 _d 0 * ( CbotCp(i,j,bi,bj)+CbotCp(i-1,j,bi,bj) )
CML     &         * utp(i,j,bi,bj) +
CML#endif
CML     &         0.5 _d 0 * ( DWATNp(i,j,bi,bj)+DWATNp(i-1,j,bi,bj) )
CML     &         * COSWAT * dup(i,j)
CML     &         - SIGN(SINWAT, _fCori(i,j,bi,bj))* 0.5 _d 0 *
CML     &         ( DWATNp(i  ,j,bi,bj) *0.5 _d 0*(dvp(i  ,j)+dvp(i  ,j+1))
CML     &         + DWATNp(i-1,j,bi,bj) *0.5 _d 0*(dvp(i-1,j)+dvp(i-1,j+1))
CML     &         ) ) * areaW(i,j)
CML          drag0 = (
CML#ifdef SEAICE_ALLOW_BOTTOMDRAG
CML     &         0.5 _d 0 * ( CbotC(i,j,bi,bj)+CbotC(i-1,j,bi,bj) )
CML     &         * uIceLoc(i,j,bi,bj) +
CML#endif
CML     &         0.5 _d 0 * ( DWATN(i,j,bi,bj)+DWATN(i-1,j,bi,bj) )
CML     &         * COSWAT * du(i,j)
CML     &         - SIGN(SINWAT, _fCori(i,j,bi,bj))* 0.5 _d 0 *
CML     &         ( DWATN(i  ,j,bi,bj) * 0.5 _d 0 * (dv(i  ,j)+dv(i  ,j+1))
CML     &         + DWATN(i-1,j,bi,bj) * 0.5 _d 0 * (dv(i-1,j)+dv(i-1,j+1))
CML     &         ) ) * areaW(i,j)
CML          uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj) + (dragp-drag0)*reps
CMLC                                      + (Cdrag*cosWat+Cb)*uIce + uIce*sinWat)
CML          dragp = (
CML#ifdef SEAICE_ALLOW_BOTTOMDRAG
CML     &         0.5 _d 0 * ( CbotCp(i,j,bi,bj)+CbotCp(i,j-1,bi,bj) )
CML     &         * vtp(i,j,bi,bj) +
CML#endif
CML     &         0.5 _d 0 * ( DWATNp(i,j,bi,bj)+DWATNp(i,j-1,bi,bj) )
CML     &         * COSWAT * dvp(i,j)
CML     &         + SIGN(SINWAT, _fCori(i,j,bi,bj)) * 0.5 _d 0 *
CML     &         ( DWATNp(i,j  ,bi,bj) *0.5 _d 0*(dup(i,j  )+dup(i+1,j  ))
CML     &         + DWATNp(i,j-1,bi,bj) *0.5 _d 0*(dup(i,j-1)+dup(i+1,j-1))
CML     &         ) ) * areaS(i,j)
CML          drag0 = (
CML#ifdef SEAICE_ALLOW_BOTTOMDRAG
CML     &         0.5 _d 0 * ( CbotC(i,j,bi,bj)+CbotC(i,j-1,bi,bj) )
CML     &         * vIceLoc(i,j,bi,bj) +
CML#endif
CML     &         0.5 _d 0 * ( DWATN(i,j,bi,bj)+DWATN(i,j-1,bi,bj) )
CML     &         * COSWAT * dv(i,j)
CML     &         + SIGN(SINWAT, _fCori(i,j,bi,bj)) * 0.5 _d 0 *
CML     &         ( DWATN(i,j  ,bi,bj) * 0.5 _d 0 * (du(i,j  )+du(i+1,j  ))
CML     &         + DWATN(i,j-1,bi,bj) * 0.5 _d 0 * (du(i,j-1)+du(i+1,j-1))
CML     &         ) ) * areaS(i,j)
CML          vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj) + (dragp-drag0)*reps
CMLC     apply masks for interior (important when we have open boundaries)
CML          uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj)*maskinW(i,j,bi,bj)
CML          vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj)*maskinS(i,j,bi,bj)
CML         ENDDO
CML        ENDDO
CML#ifdef SEAICE_ALLOW_MOM_ADVECTION
CML        IF ( SEAICEmomAdvection ) THEN
CML         DO j=1-OLy,sNy+OLy
CML          DO i=1-OLx,sNx+OLx
CML           gUmom(i,j) = 0. _d 0
CML           gVmom(i,j) = 0. _d 0
CML          ENDDO
CML         ENDDO
CML         CALL SEAICE_MOM_ADVECTION(
CML     I        bi,bj,1,sNx,1,sNy,
CML     I        utp, vtp,
CML     O        gUmom, gVmom,
CML     I        myTime, myIter, myThid )
CMLC     Beware of sign! gU/Vmom is computed for the rhs of the equation;
CMLC     therefore, we need to substract gU/Vmom from the left hand side
CML         DO j=1,sNy
CML          DO i=1,sNx
CML           uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj) - gUmom(i,j)*reps
CML           vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj) - gVmom(i,j)*reps
CML          ENDDO
CML         ENDDO
CML         CALL SEAICE_MOM_ADVECTION(
CML     I        bi,bj,1,sNx,1,sNy,
CML     I        uIceLoc, vIceLoc,
CML     O        gUmom, gVmom,
CML     I        myTime, myIter, myThid )
CML         DO j=1,sNy
CML          DO i=1,sNx
CML           uIceLHS(i,j,bi,bj) = uIceLHS(i,j,bi,bj) + gUmom(i,j)*reps
CML           vIceLHS(i,j,bi,bj) = vIceLHS(i,j,bi,bj) + gVmom(i,j)*reps
CML          ENDDO
CML         ENDDO
CML        ENDIF
CML#endif /* SEAICE_ALLOW_MOM_ADVECTION */
CML        DO j=1,sNy
CML         DO i=1,sNx
CML           duIce(i,j,bi,bj) = uIceLHS(i,j,bi,bj)
CML           dvIce(i,j,bi,bj) = vIceLHS(i,j,bi,bj)
CML         ENDDO
CML        ENDDO
CML       ENDDO
CML      ENDDO

#endif /* SEAICE_ALLOW_JFNK */

      RETURN
      END
