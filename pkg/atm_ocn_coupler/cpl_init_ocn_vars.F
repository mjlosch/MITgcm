C $Header: /u/gcmpack/MITgcm/pkg/atm_ocn_coupler/cpl_init_ocn_vars.F,v 1.1 2015/10/26 13:11:32 jmc Exp $
C $Name:  $

#include "CPP_OPTIONS.h"

CBOP 0
C !ROUTINE: CPL_INIT_OCN_VARS

C !INTERFACE:
      SUBROUTINE CPL_INIT_OCN_VARS( msgUnit )

C !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE CPL_INIT_OCN_VARS
C     | o Initialise to zero all OCNVARS.h arrays
C     *==========================================================*

C !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "OCNSIZE.h"
#include "OCNVARS.h"
c#include "OCNIDS.h"

C !INPUT PARAMETERS:
C     msgUnit          :: log-file I/O unit
      INTEGER msgUnit

C !LOCAL VARIABLES:
C     i,j     :: loop index
      INTEGER i, j
CEOP

      WRITE(msgUnit,'(2A)') 'CPL_INIT_OCN_VARS: ',
     &                      'initialise all OCNVARS.h arrays'

      DO j=1,Ny_ocn
        DO i=1,Nx_ocn
          landMask_ocn   (i,j) = 0. _d 0
          OcMxlD_ocn     (i,j) = 0. _d 0
          SST_ocn        (i,j) = 0. _d 0
          SSS_ocn        (i,j) = 0. _d 0
          vSq_ocn        (i,j) = 0. _d 0
          fluxCO2_ocn    (i,j) = 0. _d 0
          atmSLPr_ocn    (i,j) = 0. _d 0
          HeatFlux_ocn   (i,j) = 0. _d 0
          qShortWave_ocn (i,j) = 0. _d 0
          TauX_ocn       (i,j) = 0. _d 0
          TauY_ocn       (i,j) = 0. _d 0
          FWFlux_ocn     (i,j) = 0. _d 0
          SaltFlx_ocn    (i,j) = 0. _d 0
          sIceMass_ocn   (i,j) = 0. _d 0
          aCO2_ocn       (i,j) = 0. _d 0
          wSpeed_ocn     (i,j) = 0. _d 0
          sIceFrac_ocn   (i,j) = 0. _d 0
        ENDDO
      ENDDO

      RETURN
      END