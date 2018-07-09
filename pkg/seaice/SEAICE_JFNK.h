CBOP
C !ROUTINE: SEAICE_JFNK.h

C !DESCRIPTION: \bv
C     *==========================================================*
C     | SEAICE_JFNK.h
C     | o Header for JFNK solver specific fields.
C     *==========================================================*
C
C \ev
CEOP

#if (defined SEAICE_ALLOW_JFNK) || (defined SEAICE_ALLOW_KRYLOV)
      COMMON/SEAICE_DAMPED_JACOBIAN/SEAICElambdaDampedJacobian
      _RL SEAICElambdaDampedJacobian
      COMMON/SEAICE_DAMPED_JACOBIAN_FIELD/
     &     etaPre,etaZPre,zetaPre,zetaZPre
      _RL zetaPre (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zetaZPre(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etaPre  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etaZPre (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     diagnostics for the JFNK and Krylov solver
      INTEGER totalNewtonIters
      INTEGER totalNewtonFails
      INTEGER totalKrylovIters
      INTEGER totalKrylovFails
      INTEGER totalJFNKtimeSteps
      COMMON /SEAICE_SOLVER_I/
     &     totalNewtonIters, totalNewtonFails,
     &     totalKrylovIters, totalKrylovFails,
     &     totalJFNKtimeSteps
      INTEGER nVec
      PARAMETER ( nVec=2*sNx*sNy )
      _RL scalarProductMetric( nVec, 1, nSx, nSy )
      COMMON /SEAICE_KRYLOV_RL/ scalarProductMetric
#endif /* SEAICE_ALLOW_JFNK or SEAICE_ALLOW_KRYLOV */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
