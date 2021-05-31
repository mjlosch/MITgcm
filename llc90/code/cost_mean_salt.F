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
      INTEGER nMsk
      PARAMETER (nMsk=4)
      _RL locfc,tmp,tmp2,volwght,saltTmp,slvl_arc_oc,slvl_arc_coa
      _RS mask_arc_oc    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,   nSx,nSy)
      _RS mask_arc_coa   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,   nSx,nSy)
      _RL areaWghtFld_oc (1-OLx:sNx+OLx,1-OLy:sNy+OLy,   nSx,nSy)
      _RL areaWghtFld_coa(1-OLx:sNx+OLx,1-OLy:sNy+OLy,   nSx,nSy)
      _RL volWghtFld     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS maskUloc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nMsk,nSx,nSy)
      _RS maskVloc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nMsk,nSx,nSy)
      _RS hfw, hfs
      _RL dummyRL(1)
      LOGICAL exst
      character*(MAX_LEN_FNAM) fname
      integer itlo,ithi
      integer jtlo,jthi
C Functions
      INTEGER ILNBLNK, il

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

      Nk = 2
C     Read masks

      exst=.FALSE.
      fname = 'mask_arctic_ocean.data'
      INQUIRE(file = fname, exist = exst)
      IF ( exst ) THEN
       CALL READ_FLD_XY_RS(fname,' ',mask_arc_oc,0,myThid)
      ELSE
       il = ILNBLNK(fname)
       WRITE(standardmessageunit,*)
     &       'S/R COST_MEAN_SALT: ',fname(1:il), ' does not exist'
      ENDIF

      fname = 'mask_Arcitc_Ocean_coast_100km.data'
      INQUIRE(file = fname, exist = exst)
      IF ( exst ) THEN
       CALL READ_FLD_XY_RS(fname,' ',mask_arc_coa,0,myThid)
      ELSE
       il = ILNBLNK(fname)
       WRITE(standardmessageunit,*)
     &       'S/R COST_MEAN_SALT: ',fname(1:il), ' does not exist'
      ENDIF

#ifdef ALLOW_MDSIO
      fname = 'maskU_arctic_sections.data'
      INQUIRE(file = fname, exist = exst)
      IF ( exst ) THEN
       CALL MDS_READ_FIELD(
     I      fname, readBinaryPrec, .FALSE., 'RS', nMsk, 1, nMsk,
     O      dummyRL, maskUloc,
     I      1, myThid )
      ELSE
       il = ILNBLNK(fname)
       WRITE(standardmessageunit,*)
     &       'S/R COST_MEAN_SALT: ',fname(1:il), ' does not exist'
      ENDIF
      fname = 'maskV_arctic_sections.data'
      INQUIRE(file = fname, exist = exst)
      IF ( exst ) THEN
       CALL MDS_READ_FIELD(
     I      fname, readBinaryPrec, .FALSE., 'RS', nMsk, 1, nMsk,
     O      dummyRL, maskVloc,
     I      1, myThid )
      ELSE
       il = ILNBLNK(fname)
       WRITE(standardmessageunit,*)
     &       'S/R COST_MEAN_SALT: ',fname(1:il), ' does not exist'
      ENDIF
#else
      cost_mean_salt.F: Oops!
#endif

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
           locfc = locfc + volWghtFld(i,j,k,bi,bj) * saltTmp
CML           print*, 'mean salinity: ',cMeanSalt(i,j,k,bi,bj)
CML           print*, 'arctic coean mask:', mask_arc_oc(i,j,bi,bj)
CML           print*, 'logistic equation:'
CML           print*, (1-1/(1+exp(-20*(cMeanSalt(i,j,k,bi,bj)-35))))
CML           print*, 'weighted volume', volwght
CML           print*, 'mean of tile:',locfc
          ENDDO
         ENDDO
        ENDDO
        objf_mean_salt(bi,bj) = locfc * tmp
       ENDDO
      ENDDO

C     compute sea level change

C     compute area weights
      tmp  = 0. _d 0
      tmp2 = 0. _d 0
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        areaWghtFld_oc(:,:,bi,bj)  = 0. _d 0
        areaWghtFld_coa(:,:,bi,bj) = 0. _d 0
        k = 1
        DO j=1,sNy
         DO i=1,sNx
          areaWghtFld_oc (i,j,bi,bj) = mask_arc_oc(i,j,bi,bj)
     &          *hFacC(i,j,k,bi,bj)*rA(i,j,bi,bj)
          areaWghtFld_coa(i,j,bi,bj) = mask_arc_coa(i,j,bi,bj)
     &          *hFacC(i,j,k,bi,bj)*rA(i,j,bi,bj)
          tmp  = tmp  + areaWghtFld_oc (i,j,bi,bj)
          tmp2 = tmp2 + areaWghtFld_coa(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      _GLOBAL_SUM_RL( tmp , myThid )
      _GLOBAL_SUM_RL( tmp2, myThid )
      IF ( tmp .GT.0. _d 0 ) tmp  = 1. _d 0 / tmp
      IF ( tmp2.GT.0. _d 0 ) tmp2 = 1. _d 0 / tmp2

C     Compute sea level cost function
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        slvl_arc_oc  = 0. _d 0
        slvl_arc_coa = 0. _d 0
        DO j=1,sNy
         DO i=1,sNx
          slvl_arc_oc  = slvl_arc_oc  + cMeanEtaTot(i,j,bi,bj)
     &         *areaWghtFld_oc(i,j,bi,bj)
          slvl_arc_coa = slvl_arc_coa + cMeanEtaTot(i,j,bi,bj)
     &         *areaWghtFld_coa(i,j,bi,bj)
         ENDDO
        ENDDO
        objf_mean_slvl_arc_oc (bi,bj) = slvl_arc_oc  * tmp
        objf_mean_slvl_arc_coa(bi,bj) = slvl_arc_coa * tmp2
       ENDDO
      ENDDO

C     compute fluxes of salt through specific sections:
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        locfc = 0. _d 0
        objf_sflux_bering (bi,bj) = locfc
        objf_sflux_barents(bi,bj) = locfc
        objf_sflux_fram   (bi,bj) = locfc
        objf_sflux_caa    (bi,bj) = locfc
        DO k=1,Nr
         DO j=1,sNy
          DO i=1,sNx
C     convert volume flux to g/kg Sverdup = 1e6 (g/kg)(m^3/s) to deal with smaller values
           hfw = hFacW(i,j,k,bi,bj)*drF(k)*dyG(i,j,bi,bj) * 1. _d -6
           hfs = hFacS(i,j,k,bi,bj)*drF(k)*dxG(i,j,bi,bj) * 1. _d -6
C     Bering Strait
           objf_sflux_bering (bi,bj) =  objf_sflux_bering (bi,bj)
     &          + cMeanSaltUVel(i,j,k,bi,bj)*maskUloc(i,j,1,bi,bj)*hfw
     &          + cMeanSaltVVel(i,j,k,bi,bj)*maskVloc(i,j,1,bi,bj)*hfs
C     Barents Sea
           objf_sflux_barents(bi,bj) =  objf_sflux_barents(bi,bj)
     &          + cMeanSaltUVel(i,j,k,bi,bj)*maskUloc(i,j,2,bi,bj)*hfw
     &          + cMeanSaltVVel(i,j,k,bi,bj)*maskVloc(i,j,2,bi,bj)*hfs
C     Fram Strait
           objf_sflux_fram   (bi,bj) =  objf_sflux_fram   (bi,bj)
     &          + cMeanSaltUVel(i,j,k,bi,bj)*maskUloc(i,j,3,bi,bj)*hfw
     &          + cMeanSaltVVel(i,j,k,bi,bj)*maskVloc(i,j,3,bi,bj)*hfs
C     Canadian Arctic Archipelago
           objf_sflux_caa    (bi,bj) =  objf_sflux_caa    (bi,bj)
     &          + cMeanSaltUVel(i,j,k,bi,bj)*maskUloc(i,j,4,bi,bj)*hfw
     &          + cMeanSaltVVel(i,j,k,bi,bj)*maskVloc(i,j,4,bi,bj)*hfs
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C     diagnostics
      CALL WRITE_FLD_XYZ_RL('maskArcticVol',' ',volWghtFld,-1,myThid)
      mask_arc_oc = maskUloc(:,:,1,:,:)
      CALL WRITE_FLD_XY_RL('maskArcticBeringU',' ', mask_arc_oc,
     &     -1, myThid)
      mask_arc_oc = maskUloc(:,:,2,:,:)
      CALL WRITE_FLD_XY_RL('maskArcticBarentsU',' ', mask_arc_oc,
     &     -1, myThid)
      mask_arc_oc = maskVloc(:,:,2,:,:)
      CALL WRITE_FLD_XY_RL('maskArcticBarentsV',' ', mask_arc_oc,
     &     -1, myThid)
      mask_arc_oc = maskUloc(:,:,3,:,:)
      CALL WRITE_FLD_XY_RL('maskArcticFramU',' ', mask_arc_oc,
     &     -1, myThid)
      mask_arc_oc = maskUloc(:,:,4,:,:)
      CALL WRITE_FLD_XY_RL('maskArcticCaaU',' ', mask_arc_oc,
     &     -1, myThid)
      mask_arc_oc = maskVloc(:,:,4,:,:)
      CALL WRITE_FLD_XY_RL('maskArcticCaaV',' ', mask_arc_oc,
     &     -1, myThid)
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        write(standardmessageunit,*)
     &       'ml-objf_mean_salt     =',objf_mean_salt(bi,bj),
     &       tmp,startTime,endTime
        write(standardmessageunit,*)
     &       'ml-objf_sflux_bering      =', objf_sflux_bering (bi,bj)
        write(standardmessageunit,*)
     &       'ml-objf_sflux_barents     =', objf_sflux_barents(bi,bj)
        write(standardmessageunit,*)
     &       'ml-objf_sflux_fram        =', objf_sflux_fram   (bi,bj)
        write(standardmessageunit,*)
     &       'ml-objf_sflux_caa         =', objf_sflux_caa    (bi,bj)
        write(standardmessageunit,*)
     &       'ml-objf_mean_slvl_arc_oc  =', objf_mean_slvl_arc_oc(bi,bj)
        write(standardmessageunit,*)
     &       'ml-objf_mean_slvl_arc_coa =',objf_mean_slvl_arc_coa(bi,bj)
       ENDDO
      ENDDO

#endif /* ALLOW_COST_MEAN_SALT */
      RETURN
      END
