#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#else
# define OBCS_UVICE_OLD
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C--  File seaice_update_damage.F: seaice LSR dynamical solver S/R:
C--   Contents
C--   o SEAICE_UPDATE_DAMAGE
C--   o SEAICE_RECIP_RELAXTIME
C--   o SEAICE_MEB_UPDATE_RHS
C--   o SEAICE_MEB_CALC_SIGMA
C--   o SEAICE_MEB_CALC_STRESS
C--   o SEAICE_MEB_UPDATE_SIGMA

CBOP
C     !ROUTINE: SEAICE_UPDATE_DAMAGE
C     !INTERFACE:
      SUBROUTINE SEAICE_UPDATE_DAMAGE(
     &     updateSigma, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_UPDATE_DAMAGE
C     | o determine new damage based on present stress field
C     |   V. Dansereau et al, 2016, TC
C     *==========================================================*
C     | written by Martin Losch, Jun 2017
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
#include "SEAICE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      LOGICAL updateSigma

#if ( defined (SEAICE_CGRID) && defined (SEAICE_ALLOW_DYNAMICS) )
# ifdef SEAICE_ALLOW_DAMAGE

C     === Local variables ===
C     i,j,bi,bj      :: Loop counters
C     e11,e12,e22    :: components of strain rate tensor
C     seaice_div     :: divergence strain rates at C-points times P
C                       /divided by Delta minus 1
C     seaice_tension :: tension    strain rates at C-points times P
C                       /divided by Delta
C     seaice_shear   :: shear      strain rates, defined at Z-points times P
C                       /divided by Delta
C     sig11, sig22   :: sum and difference of diagonal terms of stress tensor
C     modification for adaptive alpha and beta
C               (see Kimmritz, Danilov, Losch 2015 for gamma << alpha beta)
C     EVPcFac        ::SEAICE_deltaTdyn*SEAICEaEVPcStar*(SEAICEaEVPcoeff*PI)**2
C                        with
C     SEAICEaEVPcStar:: multiple of stabilty factor: alpha*beta = cstar*gamma
C     SEAICEaEvpcoeff:: largest stabilized frequency according to
C                        gamma = zeta * (cfac/cellarea)*deltaT/m
C                                with   (cfac/cellarea) <= pi**2/cellarea
C     evpAlphaC/Z    :: alpha field on C points and on Z points
C                        := sqrt(cstar gamma)
C     evpBetaU/V     :: beta field on u and on v points
C                        := sqrt(cstar gamma)
C     evpAlphaMin    :: lower limit of alpha and beta, regularisation
C                     to prevent singularities of system matrix,
C                     e.g. when ice concentration is too low.
C     betaFacP1U/V   :: = betaFacP1 in standard case,
C                          with varying beta in the adaptive case
C                           on u and on v point
C     betaFacU/V     :: analog betaFacP1U/V

      INTEGER i, j, bi, bj

      _RL deltaTloc
      _RL q, recip_q, sqrtq
      _RL damageCrit
      _RL sumNorm
      _RL recip_healingTime, recip_damageTime
      _RL sig1, sig2, sigTmp, sigp, sigm
      _RL dtmp(1-OLx:sNx+OLx,1-OLy:sNy+OLy), wght

      _RL sigmaT        (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaC        (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig11         (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig22         (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12         (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12C        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL MC_FAILURE, TENSILE_FAILURE, fldAtZ, fldAtC
      EXTERNAL MC_FAILURE, TENSILE_FAILURE, fldAtZ, fldAtC

CEOP

C
      deltaTloc = SEAICE_deltaTdyn

C     parameters
      recip_healingTime= 0.
      IF ( SEAICEhealingTime .GT. 0. _d 0 )
     &     recip_healingTime = 1. _d 0/SEAICEhealingTime
      recip_damageTime = 0.
      IF ( SEAICEdamageTime .GT. 0. _d 0 )
     &     recip_damageTime  = 1. _d 0/SEAICEdamageTime
C     slope of the Mohr-Coulomb-envelope in principle stress plane
      sqrtq = SQRT(SEAICEintFrictCoeff**2 + 1) + SEAICEintFrictCoeff
C     sqrtq = SQRT(0.5*SEAICEmohrCoulombSlope/(1.-SEAICEmohrCoulombSlope**2)
      q = sqrtq**2
      recip_q = 1. _d 0 / q
C     compute latest strain rates
      CALL SEAICE_CALC_STRAINRATES(
     I     uIce, vIce,
     O     e11, e22, e12,
     I     0, myTime, myIter, myThid )
C
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO J=1-OLy,sNy+OLy
         DO I=1-OLx,sNx+OLx
C     uniaxial (unconfined) compressive strength (needs to include HEFF
C     because of our definition of stress including HEFF) is negative
C     according to our convention (this is different from Dansereau et al 2016)
          sigmaC(I,J) = - 2.*cohesion(I,J,bi,bj)*sqrtq*HEFF(I,J,bi,bj)
C     tensile strength cutoff
          sigmaT(I,J) = - sigmaC(I,J)*recip_q
C
          dCrit(I,J,bi,bj)  = 1. _d 0*maskInC(I,J,bi,bj)
          dCritZ(I,J,bi,bj) = 1. _d 0
         ENDDO
        ENDDO

C     re-compute stress from most recent strain rates
#  ifdef SEAICE_ALLOW_MEB
        IF ( SEAICEuseMEB ) THEN
         CALL SEAICE_MEB_CALC_STRESS(
     I        e11, e22, e12, press0, etaZ, damage,
     I        seaice_sigma1, seaice_sigma2, seaice_sigma12,
     I        seaice_sigma12C,
     O        sig11, sig22, sig12, sig12C,
     I        bi, bj, myTime, myIter, myThid )
        ELSE
#  else
        IF ( .TRUE. ) THEN
#  endif /* not SEAICE_ALLOW_MEB */
         CALL SEAICE_CALC_STRESS(
     I        e11, e22, e12, press, zeta, eta, etaZ,
     O        sig11, sig22, sig12,
     I        bi, bj, myTime, myIter, myThid )
        ENDIF
        DO j=1,sNy
         DO i=1,sNx
C     abbreviations
          sigp   = sig11(I,J) + sig22(I,J)
          sigm   = sig11(I,J) - sig22(I,J)
          sigTmp = SQRT( sigm*sigm + 4.*sig12C(I,J,bi,bj)**2 )
C     principle components
          sig1   = 0.5 _d 0 * ( sigp + sigTmp )
          sig2   = 0.5 _d 0 * ( sigp - sigTmp )
C     damageCrit = 1 means not additional damage
          damageCrit = 1. _d 0
C     here we can clearly do something about regularisation
C     According to our convention, tensile stress is positive and our
C     Mohr-Coulomb criterion reads sig2-sig1*q = sigmaC with negative sigmaC.
C     The tensile failure criterion switch sig1 and sig2 with respect to
C     Dansereau et al (2016) so that for sig1>0, sig2>sigmaT
CML          IF ( SEAICEuseTensileCutoff .AND.
CML     &         sig1 .GT. 0. _d 0 .AND. sig2 .GT. sigmaT(I,J) )
          damageCrit = TENSILE_FAILURE(sig1,sig2,sigmaT(I,J),damageCrit)
          damageCrit = TENSILE_FAILURE(sig2,sig1,sigmaT(I,J),damageCrit)
C     keep in mind that sigmaC < 0 so that the following criterion catches
C     the case when sig2 is below sigmaC+q*sig1
          damageCrit = MC_FAILURE(sig1,sig2,q,sigmaC(I,J),damageCrit)
C     or above sigmaC+q*sig2
          damageCrit = MC_FAILURE(sig2,sig1,q,sigmaC(I,J),damageCrit)
C     critical damage cannot be below zero (eq 39 in mebdoc.tex)
          damageCrit = MAX(damageCrit, 0. _d 0)
C Damage-Options:
C Option A1: Damage    explicit and for dt< damagetime:
          damage(I,J,bi,bj) = damageNm1(I,J,bi,bj) *
     &         (deltaTloc*recip_damageTime*(damageCrit-1.)+ 1.)
     &          + deltaTloc*recip_healingTime * maskInC(I,J,bi,bj)
C Option A2: new old version of damage without dependeny of time step
c          damage(I,J,bi,bj) = damageNm1(I,J,bi,bj)*damageCrit(I,J,bi,bj) 
c          damage(I,J,bi,bj) = damage(I,J,bi,bj) +
c     &     (1-damage(I,J,bi,bj)*deltaTLoc*recip_healingTime)
C Option A3: do not update damage ever
c          damage(i,j,bi,bj) = 1.
C     make sure that 0 < damage <=1
          damage(I,J,bi,bj) = MIN(1. _d 0, damage(I,J,bi,bj))
          damage(I,J,bi,bj) = MAX(0. _d 0, damage(I,J,bi,bj))
C Sigma-Update-Options:
C Option B1: update stress as time dependent as dynamic parameters
C          dCrit(I,J,bi,bj) = damage(I,J,bi,bj)
C Option B2: update with dcrit
          dCrit(I,J,bi,bj) = damageCrit*maskInC(I,J,bi,bj)
C Option B3: never update dcrit
c          dCrit(I,J,bi,bj) = maskInC(I,J,bi,bj)
          IF ( AREA(I,J,bi,bj) .EQ. 0. _d 0 )
     &         damage(I,J,bi,bj) = 0. _d 0
          IF ( updateSigma ) THEN
C     compute updated/reduced stress that lies on or within yield curve
          sig12C(I,J,bi,bj)  = sig12C(I,J,bi,bj) * dCrit(I,J,bi,bj)
          seaice_sigma1(I,J,bi,bj) = sigp        * dCrit(I,J,bi,bj)
          seaice_sigma2(I,J,bi,bj) = sigm        * dCrit(I,J,bi,bj)
          seaice_sigma12C(I,J,bi,bj)=sig12C(I,J,bi,bj)
          ENDIF
C     just diagnostics
C     hack to save principle stress before updating sigma
          sigTmp = 0.
          IF ( HEFF(I,J,bi,bj).GT.1. _d -4)
     &         sigTmp=1./SEAICE_strength/HEFF(I,J,bi,bj)
          sig11(I,J)  = sig1*sigTmp
          sig22(I,J)  = sig2*sigTmp
         ENDDO
        ENDDO
        CALL DIAGNOSTICS_FILL(sig11,'SIsig1pr',0,1,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(sig22,'SIsig2pr',0,1,2,bi,bj,myThid)
C dcritZ-Option
C Option C1/C2: averaged/new dcritz
        DO j=1,sNy+1
         DO i=1,sNx+1
          sigp   = fldAtZ(sig11 + sig22,i,j,bi,bj)
          sigm   = fldAtZ(sig11 - sig22,i,j,bi,bj)
          sigTmp = SQRT( sigm*sigm + 4.*sig12(i,j)**2 )
          sig1   = 0.5 _d 0 * ( sigp + sigTmp )
          sig2   = 0.5 _d 0 * ( sigp - sigTmp )
          damageCrit = 1. _d 0
          sigmaC(I,J) = - 2.*cohesion(I,J,bi,bj)*sqrtq
     &         *fldAtZ(HEFF(:,:,bi,bj),i,j,bi,bj)
          sigmaT(I,J) = - sigmaC(I,J)*recip_q
          damageCrit = TENSILE_FAILURE(sig1,sig2,sigmaT(I,J),damageCrit)
          damageCrit = TENSILE_FAILURE(sig2,sig1,sigmaT(I,J),damageCrit)
          damageCrit = MC_FAILURE(sig1,sig2,q,sigmaC(I,J),damageCrit)
          damageCrit = MC_FAILURE(sig2,sig1,q,sigmaC(I,J),damageCrit)
          damageCrit = MAX(damageCrit, 0. _d 0)
C Damage-Options:
C Option A1: Damage    explicit and for dt< damagetime:
          damageZ(I,J,bi,bj) = damageNm1Z(I,J,bi,bj) *
     &         (deltaTloc*recip_damageTime*(damageCrit-1.)+ 1.)
     &          + deltaTloc*recip_healingTime
C Option A2: damage without dependency of time step
c          damageZ(I,J,bi,bj) = damageNm1Z(I,J,bi,bj)*damageCrit(I,J,bi,bj) 
c          damageZ(I,J,bi,bj) = damageZ(I,J,bi,bj) +
c     &     (1-damage(I,J,bi,bj)*deltaTLoc*recip_healingTime)
C Option A3: do not update damage ever
c          damageZ(i,j,bi,bj) = 1.
C     make sure that 0 < damage <=1
          damageZ(I,J,bi,bj) = MIN(1. _d 0, damageZ(I,J,bi,bj))
          damageZ(I,J,bi,bj) = MAX(0. _d 0, damageNm1Z(I,J,bi,bj))
C Sigma-Update-Options:
C Option B1/C1: update stress time dependent/averaged 
C          dCritZ(I,J,bi,bj) = fldAtZ(damageZ(:,:,bi,bj),I,J,bi,bj)
C Option B1/C2: update stress time dependent/new 
C          dCritZ(I,J,bi,bj) = damageZ(I,J,bi,bj) 
C Option B2/C1: update with dcrit/averaged
C          dCritZ(i,j,bi,bj) = fldAtZ(dCrit(:,:,bi,bj),I,J,bi,bj)
C Option B2/C2: update with new dcrit
          dCritZ(I,J,bi,bj) = damageCrit
C option B3: never update dcrit
c      dCritZ(I,J,bi,bj) = 1. _d 0
         ENDDO
        ENDDO  
C If out of Loop for better performance
        IF ( updateSigma ) THEN
         DO J=1, sNy+1
          DO I=1, sNx+1
           seaice_sigma12(I,J,bi,bj) = sig12(I,J)
     &         * dCritZ(I,J,bi,bj)
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
      _EXCH_XY_RL(sig12C, myThid )
      _EXCH_XY_RL(dCrit, myThid )
      _EXCH_XY_RL(damage, myThid )
C
# endif /* SEAICE_ALLOW_DAMAGE */

      RETURN
      END

      _RL FUNCTION MC_FAILURE ( sig1, sig2, q, sigC, damageCrit )
      IMPLICIT NONE
C     compute Mohr-Coulomb failure criterion and updage damageCrit
      _RL sig1, sig2, q, sigC, damageCrit
      _RL denom, rdenm
      denom = sig2 - q*sig1
      rdenm = 0. _d 0
      IF ( denom .NE. 0. _d 0 ) rdenm = 1./denom
      mc_failure = damageCrit
      IF ( denom .LT. sigC ) mc_failure = MIN(damageCrit,sigC*rdenm)
      RETURN
      END
C
      _RL FUNCTION TENSILE_FAILURE ( sig1, sig2, sigT, damageCrit )
      IMPLICIT NONE
C     compute tensile failure criterion and updage damageCrit
      _RL sig1, sig2, sigT, damageCrit
      tensile_failure = damageCrit
      IF ( sig1 .GT. 0. _d 0 .AND. sig2 .GT. sigT )
     &     tensile_failure = MIN(damageCrit,sigT/sig2)
      RETURN
      END
C Averaging-Options F1/F2/F3/F4
      _RL FUNCTION fldAtZ ( fldIn, i, j, bi, bj )
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
      INTEGER i,j, bi,bj
      _RL fldIn (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sumNorm
C
      sumNorm  = maskInC(I,J,  bi,bj)+maskInC(I-1,J,  bi,bj)
     &     +     maskInC(I,J-1,bi,bj)+maskInC(I-1,J-1,bi,bj)
      IF ( sumNorm.GT.0. _d 0 ) sumNorm = 1. _d 0 / sumNorm
C Option-G1: simple mean on C-Points
C      fldAtZ = sumNorm *
C     &     ( fldIn(I,  J  )*maskInC(I,  J,  bi,bj)
C     &     + fldIn(I-1,J  )*maskInC(I-1,J,  bi,bj)
C     &     + fldIn(I,  J-1)*maskInC(I  ,J-1,bi,bj)
C     &     + fldIn(I-1,J-1)*maskInC(I-1,J-1,bi,bj) )
C Option-G2: quadratic mean on C-Points
      fldAtZ = SQRT (sumNorm *
     &     ( (fldIn(I,  J  )*maskInC(I,  J,  bi,bj))**2
     &     + (fldIn(I-1,J  )*maskInC(I-1,J,  bi,bj))**2
     &     + (fldIn(I,  J-1)*maskInC(I  ,J-1,bi,bj))**2
     &     + (fldIn(I-1,J-1)*maskInC(I-1,J-1,bi,bj))**2 ))
C
      RETURN
      END

      _RL FUNCTION fldAtC ( fldIn, i, j, bi,bj )
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
      INTEGER i,j, bi,bj
      _RL fldIn (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sumNorm
C
      sumNorm  = 0.25 _d 0
C Option-G3: simple mean on Z-Points
      fldAtC = (sumNorm * recip_rA(I,J,bi,bj) *
     &     ( rAz(I,  J,  bi,bj)*fldIn(I,  J  )
     &     + rAz(I+1,J,  bi,bj)*fldIn(I+1,J  )
     &     + rAz(I,  J+1,bi,bj)*fldIn(I,  J+1)
     &     + rAz(I+1,J+1,bi,bj)*fldIn(I+1,J+1)
     &     ) ) * maskInC(I,J,bi,bj)
C Option-G4: quadratic mean on C-Points
C      fldAtC = SQRT( sumNorm * recip_rA(I,J,bi,bj) *
C     &     ( rAz(I,  J,  bi,bj)*fldIn(I,  J  )**2
C     &     + rAz(I+1,J,  bi,bj)*fldIn(I+1,J  )**2
C     &     + rAz(I,  J+1,bi,bj)*fldIn(I,  J+1)**2
C     &     + rAz(I+1,J+1,bi,bj)*fldIn(I+1,J+1)**2
C     &     ) ) * maskInC(I,J,bi,bj)
C
      RETURN
      END

CBOP
C     !ROUTINE: SEAICE_RECIP_RELAXTIME
C     !INTERFACE:
      SUBROUTINE SEAICE_RECIP_RELAXTIME(
     I     damage,
     O     recip_rlx,
     I     bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_UPDATE_DAMAGE
C     | o determine relaxation time scale based on new damage
C     |   V. Dansereau et al, 2016, TC
C     *==========================================================*
C     | written by Martin Losch, Jun 2017
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
      _RL damage   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL recip_rlx(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

# ifdef SEAICE_ALLOW_MEB
C     === Local variables ===
C     recip_relaxTime :: ratio of undamaged elasticity modululs and
C                        undamaged appararent viscosity
      INTEGER i, j
      _RL     recip_relaxTime
CEOP

C     zero means zero viscosity, hence elasto-brittle wihout the
C     Maxwell term
      recip_relaxTime = 0. _d 0
      IF ( SEAICEviscosity .GT. 0. _d 0 )
     &     recip_relaxTime = SEAICE_strength/SEAICEviscosity
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        recip_rlx(i,j) = 0. _d 0
C     faster for damage**(1-SEAICEdamageParm), but requires damage > 0
C     This is probably also a place for regularisation?
CML        IF ( damage(i,j,bi,bj) .GT. SEAICEdamageMin )
CML     &       recip_rlx(i,j) = recip_relaxTime
CML     &       * EXP((1.-SEAICEdamageParm)*LOG(damage(i,j,bi,bj)))
        IF ( damage(i,j,bi,bj) .GT. SEAICEdamageMin )
     &       recip_rlx(i,j) = recip_relaxTime
     &       *damage(i,j,bi,bj)**(1.-SEAICEdamageParm)
        damage(i,j,bi,bj) = damage(i,j,bi,bj)*maskInC(i,j,bi,bj)
        recip_rlx(i,j)    = recip_rlx(i,j)   *maskInC(i,j,bi,bj)
       ENDDO
      ENDDO

      RETURN
      END

CBOP
C     !ROUTINE: SEAICE_MEB_UPDATE_RHS
C     !INTERFACE:
      SUBROUTINE SEAICE_MEB_UPDATE_RHS(
     I     damage, sigPlus, sigMinus, sigma12,
     U     uIceRHS, vIceRHS,
     I     bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_MEB_UPDATE_RHS
C     | o update right hand side with stress divergence of the
C     |   previous time step
C     *==========================================================*
C     | written by Martin Losch, Jun 2017
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
      _RL damage  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigPlus (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigMinus(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigma12 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     u/vIceRHS :: RHS of momentum equations
      _RL uIceRHS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceRHS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     === Local variables ===
C     i,j :: loop indices
      INTEGER i,j
      _RL sig11(0:sNx,0:sNy)
      _RL sig22(0:sNx,0:sNy)
      _RL sig12(1:sNx+1,1:sNy+1)
      _RL recip_rlx(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL denom(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sumNorm, deltaTloc
      _RL fldAtZ
      EXTERNAL fldAtZ
CEOP

C     add stress divergence of previous timestep
      deltaTloc = SEAICE_deltaTdyn
      CALL SEAICE_RECIP_RELAXTIME(
     I     damage,
     O     recip_rlx,
     I     bi, bj, myTime, myIter, myThid )
C     recover sigma11 and sigma22 and scale by denominator
C     (1 + dt*Elasiticy/viscosity)
      DO J=1-OLy,sNy+OLy
       DO I=1-OLx,sNx+OLx
        denom(I,J) = 1. _d 0 /(1. _d 0+deltaTloc*recip_rlx(I,J))
       ENDDO
      ENDDO
      DO j=0,sNy
       DO i=0,sNx
        sig11(I,J) = 0.5 _d 0 * denom(I,J) *
     &       ( sigPlus(I,J,bi,bj)+sigMinus(I,J,bi,bj) )
        sig22(I,J) = 0.5 _d 0 * denom(I,J) *
     &       ( sigPlus(I,J,bi,bj)-sigMinus(I,J,bi,bj) )
       ENDDO
      ENDDO
      DO J=1,sNy+1
       DO I=1,sNx+1
C Update-RHS-Options E1/E2
C Option E1: average whole denominator
        sig12(I,J) = sigma12(I,J,bi,bj)*fldAtZ(denom,i,j,bi,bj)
C Option E2: average only recip_rlx
C        sig12(I,J) = sigma12(I,J,bi,bj)
C     &       /(1. _d 0+deltaTloc*fldAtZ(recip_rlx,i,j,bi,bj))
       ENDDO
      ENDDO
      IF (.NOT.SEAICE_no_slip) THEN
       DO J=1,sNy+1
        DO I=1,sNx+1
         sig12(I,J) =  sig12(I,J)
     &        *maskInC(I,J,  bi,bj)*maskInC(I-1,J,  bi,bj)
     &        *maskInC(I,J-1,bi,bj)*maskInC(I-1,J-1,bi,bj)
        ENDDO
       ENDDO
      ENDIF
      DO J=1,sNy
       DO I=1,sNx
        uIceRHS(I,J,bi,bj) = uIceRHS(I,J,bi,bj) +
     &       ( sig11(I  ,J) * _dyF(I  ,J,bi,bj)
     &       - sig11(I-1,J) * _dyF(I-1,J,bi,bj)
     &       + sig12(I,J+1) * _dxV(I,J+1,bi,bj)
     &       - sig12(I,J  ) * _dxV(I,J  ,bi,bj)
     &       ) * recip_rAw(I,J,bi,bj)
        vIceRHS(I,J,bi,bj) = vIceRHS(I,J,bi,bj) +
     &       ( sig22(I,J  ) * _dxF(I,J  ,bi,bj)
     &       - sig22(I,J-1) * _dxF(I,J-1,bi,bj)
     &       + sig12(I+1,J) * _dyU(I+1,J,bi,bj)
     &       - sig12(I  ,J) * _dyU(I  ,J,bi,bj)
     &       ) * recip_rAs(I,J,bi,bj)
       ENDDO
      ENDDO

      RETURN
      END
CBOP
C     !ROUTINE: SEAICE_MEB_CALC_SIGMAP
C     !INTERFACE:
      SUBROUTINE SEAICE_MEB_CALC_SIGMAP(
     I     e11, e22, e12, press0, etaZ, damage, dCrit,
     U     sigPlus, sigMinus, sigma12, sigma12C,
     I     bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_MEB_CALC_SIGMAP
C     | o update global variables sigPlus, sigMinus, sigma12
C     *==========================================================*
C     | written by Martin Losch, Jun 2017
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
      _RL damage  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dCrit   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigPlus (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigMinus(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigma12 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CMB Plante (2020)
      _RL sigma12C(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e11     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etaZ    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     === Local variables ===
C     i,j :: loop indices
      INTEGER i,j
      _RL sig11   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig22   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12C  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CEOP

C     update press0 to include new damage
      CALL SEAICE_CALC_ICE_STRENGTH(
     I     bi, bj, myTime, myIter, myThid )
C
      CALL SEAICE_MEB_CALC_STRESS(
     I     e11, e22, e12, press0, etaZ, damage,
     I     sigPlus, sigMinus, sigma12, sigma12C,
     O     sig11, sig22, sig12, sig12C,
     I     bi, bj, myTime, myIter, myThid )

      DO J=0,sNy
       DO I=0,sNx
        sigPlus (I,J,bi,bj) = sig11(I,J) + sig22(I,J)
        sigMinus(I,J,bi,bj) = sig11(I,J) - sig22(I,J)
        sigma12C(I,J,bi,bj) = sig12C(I,J)
       ENDDO
      ENDDO
      DO J=1,sNy+1
       DO I=1,sNx+1
        sigma12(I,J,bi,bj) = sig12(I,J)
       ENDDO
      ENDDO
      RETURN
      END

CBOP
C     !ROUTINE: SEAICE_MEB_CALC_STRESS
C     !INTERFACE:
      SUBROUTINE SEAICE_MEB_CALC_STRESS(
     I     e11, e22, e12, press0, etaZ, damage,
     I     sigPlus, sigMinus, sigma12, sigma12C,
     O     sig11, sig22, sig12, sig12C,
     I     bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_MEB_CALC_STRESS
C     | o calculate new stress field from current ice state
C     |   with up-to-date damage
C     *==========================================================*
C     | written by Martin Losch, Jun 2017
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
      _RL damage  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigPlus (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigMinus(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigma12 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigma12C (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e11     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etaZ    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     ouptut :: new stress tensor components
      _RL sig11   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig22   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12C  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     === Local variables ===
C     i,j :: loop indices
      INTEGER i,j
      _RL recip_rlx(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL denom(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL denom_z
      _RL sumNorm, deltaTloc
      _RL SEAICEmebMu, SEAICEmebLambda
      _RL zetaLoc, etaLoc
      _RL sigp, sigm, ep, em, rhs12
      _RL fldAtZ, fldAtC
      EXTERNAL fldAtZ, fldAtC
CEOP

      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        sig11(i,j) = 0. _d 0
        sig22(i,j) = 0. _d 0
        sig12(i,j) = 0. _d 0
        sig12C(i,j,bi,bj) = 0. _d 0
       ENDDO
      ENDDO
C     some abbreviations
      deltaTloc = SEAICE_deltaTdyn
C     derived parameters (should probably be global variables)
      SEAICEmebMu = 0.5 _d 0/(1.+SEAICEpoissonRatio)
      SEAICEmebLambda = SEAICEmebMu * 2. _d 0*SEAICEpoissonRatio
     &     /(1. _d 0 - SEAICEpoissonRatio)
C
      CALL SEAICE_RECIP_RELAXTIME(
     I     damage,
     O     recip_rlx,
     I     bi, bj, myTime, myIter, myThid )
C
C     step stress equations
C
      DO J=1-OLy,sNy+OLy
       DO I=1-OLx,sNx+OLx
        denom(I,J) = 1. _d 0 /(1. _d 0+deltaTloc*recip_rlx(I,J))
       ENDDO
      ENDDO
      DO J=1-OLy,sNy+OLy-1
       DO I=1-OLx,sNx+OLx-1
        ep   = e11(i,j,bi,bj) + e22(i,j,bi,bj)
        em   = e11(i,j,bi,bj) - e22(i,j,bi,bj)
        zetaLoc = (SEAICEmebLambda + SEAICEmebMu) * press0(I,J,bi,bj)
        etaLoc  = SEAICEmebMu * press0(I,J,bi,bj)
        sigp    = maskInC(I,J,bi,bj)
     &       * (sigPlus (I,J,bi,bj)+deltaTloc*2.*zetaLoc*ep)*denom(I,J)
        sigm    = maskInC(I,J,bi,bj)
     &       * (sigMinus(I,J,bi,bj)+deltaTloc*2.* etaLoc*em)*denom(I,J)
C     recover sigma11 and sigma22
        sig11(I,J) = 0.5 _d 0 * ( sigp+sigm )
        sig22(I,J) = 0.5 _d 0 * ( sigp-sigm )
CMB   averaging
        rhs12 = 2.*etaLoc*fldAtC(e12(:,:,bi,bj),i,j,bi,bj)
C Calc-Stress-Options D1/D2
C Option D1: averaging sigma12
        sig12C(I,J,bi,bj) = denom(I,J) * ( deltaTLoc*rhs12 +
     &       fldAtC(sigma12(:,:,bi,bj),i,j,bi,bj) )
C Option D2: using uncoupled sigma12C
c        sig12C(I,J,bi,bj) = denom(I,J)
c     &   *(deltaTLoc*rhs12+sigma12C(I,J,bi,bj))
       ENDDO
      ENDDO
      DO J=1-OLy+1,sNy+OLy
            DO I=1-OLx+1,sNx+OLx
C     averaging to Z-points
        etaLoc = SEAICEmebMu *fldAtZ(press0(:,:,bi,bj),i,j,bi,bj)
        rhs12 = 2. _d 0 * etaLoc * e12(I,J,bi,bj)
C denom2-Options: F1/F2
C Option F1: averaging whole denominator
        denom_z = fldAtZ(denom,i,j,bi,bj)
C Option F2: averaging only recip_rlx
CML        denom_z = 1. _d 0
CML     &       /(1. _d 0+deltaTloc*fldAtZ(recip_rlx,i,j,bi,bj))
        sig12(I,J) = (sigma12(I,J,bi,bj)+deltaTLoc*rhs12)*denom_z
       ENDDO
      ENDDO
      IF (.NOT.SEAICE_no_slip) THEN
       DO J=1-OLy+1,sNy+OLy
        DO I=1-OLx+1,sNx+OLx
         sig12(I,J) =  sig12(I,J)
     &        *maskInC(I,J,  bi,bj)*maskInC(I-1,J,  bi,bj)
     &        *maskInC(I,J-1,bi,bj)*maskInC(I-1,J-1,bi,bj)
        ENDDO
       ENDDO
      ENDIF

# endif /* SEAICE_ALLOW_MEB */
#endif /* SEAICE_ALLOW_DYNAMICS and SEAICE_CGRID */
      RETURN
      END
CBOP
C     !ROUTINE: SEAICE_MEB_UPDATE_SIGMA
C     !INTERFACE:
      SUBROUTINE SEAICE_MEB_UPDATE_SIGMA(
     I     bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_MEB_UPDATE_SIGMA
C     | o update global variables sigPlus, sigMinus, sigma12 with
C     |   damage
C     *==========================================================*
C     | written by Mirjam Bouregtt, Feb 2022
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
#include "SEAICE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
C
CMB      _RL dCrit    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigPlus  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigMinus (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigma12  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sigma12C (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     === Local variables ===
C     i,j :: loop indices
      INTEGER I,J
      _RL fldAtZ, fldAtC
      _RL wght
      EXTERNAL fldAtZ, fldAtC
CEOP
      DO J=0,sNy
       DO I=0,sNx
        seaice_sigma1(I,J,bi,bj)  =seaice_sigma1(I,J,bi,bj)
     &         * dCrit(I,J,bi,bj)
        seaice_sigma2(I,J,bi,bj)  =seaice_sigma2(I,J,bi,bj)
     &         * dCrit(I,J,bi,bj)
C option D2: uncoupled sig12C
c        seaice_sigma12C(I,J,bi,bj)=seaice_sigma12C(I,J,bi,bj)
c     &         * dCrit(I,J,bi,bj)
       ENDDO
      ENDDO
      DO J=1,sNy+1
       DO I=1,sNx+1
        seaice_sigma12(I,J,bi,bj) =seaice_sigma12(I,J,bi,bj)
     &       * dCritZ(i,j,bi,bj)
       ENDDO
      ENDDO
C      wght=1. _d 0
C      DO J=1,sNy
C       DO I=1,sNx
C        seaice_sigma12C(I,J,bi,bj) =
C     &       (1. _d 0 - wght)*seaice_sigma12C(I,J,bi,bj)
C     &       + wght * fldAtC(seaice_sigma12(:,:,bi,bj),I,J,bi,bj)
C       ENDDO
C      ENDDO

      RETURN
      END
