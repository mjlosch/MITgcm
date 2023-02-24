CBOP
C !ROUTINE: CTRL_OPTIONS.h
C !INTERFACE:
C #include "CTRL_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for Control (ctrl) package:
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifndef CTRL_OPTIONS_H
#define CTRL_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_CTRL
#ifdef ECCO_CPPOPTIONS_H

C-- When multi-package option-file ECCO_CPPOPTIONS.h is used (directly included
C    in CPP_OPTIONS.h), this option file is left empty since all options that
C   are specific to this package are assumed to be set in ECCO_CPPOPTIONS.h

#else /* ndef ECCO_CPPOPTIONS_H */
C   ==================================================================
C-- Package-specific Options & Macros go here

#undef EXCLUDE_CTRL_PACK
#define ALLOW_NONDIMENSIONAL_CONTROL_IO
#define CTRL_SET_PREC_32
#define ALLOW_PACKUNPACK_METHOD2

C       >>> Other Control.
#define ALLOW_DIFFKR_CONTROL
#define ALLOW_KAPGM_CONTROL
#define ALLOW_KAPREDI_CONTROL
#define ALLOW_BOTTOMDRAG_CONTROL

C       >>> Generic Control.
#define ALLOW_GENARR2D_CONTROL
#define ALLOW_GENARR3D_CONTROL
#define ALLOW_GENTIM2D_CONTROL

C  o Rotation of wind/stress controls adjustments
C    from Eastward/Northward to model grid directions
#define ALLOW_ROTATE_UV_CONTROLS

C  o Originally the first two time-reccords of control
C    variable tau u and tau v were skipped.
C    The CTRL_SKIP_FIRST_TWO_ATM_REC_ALL option extends this
C    to the other the time variable atmospheric controls.
#undef CTRL_SKIP_FIRST_TWO_ATM_REC_ALL

C  Note: this flag turns on extra smoothing code in ctrl_get_gen.F which
C  is inconsistent with the Weaver and Courtier, 2001 algorithm, and
C  should probably not be used. The corresponding 3D flag applied only
C  to deprecated code that is now removed. At some point we will remove
C  this flag and associated code as well.
C  o apply pkg/smooth/smooth_diff2d.F to 2D controls (outside of ctrlSmoothCorrel2D)
#undef ALLOW_SMOOTH_CTRL2D

C  o use pkg/smooth correlation operator (incl. smoother) for 3D controls (Weaver, Courtier 01)
C    This CPP option just sets the default for ctrlSmoothCorrel23 to .TRUE.
c#define ALLOW_SMOOTH_CORREL3D
C  o use pkg/smooth correlation operator (incl. smoother) for 2D controls (Weaver, Courtier 01)
C    This CPP option just sets the default for ctrlSmoothCorrel2D to .TRUE.
c#define ALLOW_SMOOTH_CORREL2D

C   ==================================================================
#endif /* ndef ECCO_CPPOPTIONS_H */
#endif /* ALLOW_CTRL */
#endif /* CTRL_OPTIONS_H */
