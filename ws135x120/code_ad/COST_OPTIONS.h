CBOP
C !ROUTINE: COST_OPTIONS.h
C !INTERFACE:
C #include "COST_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for Cost-Function (cost) package:
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifndef COST_OPTIONS_H
#define COST_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_COST
#ifdef ECCO_CPPOPTIONS_H

C-- When multi-package option-file ECCO_CPPOPTIONS.h is used (directly included
C    in CPP_OPTIONS.h), this option file is left empty since all options that
C   are specific to this package are assumed to be set in ECCO_CPPOPTIONS.h

#else /* ndef ECCO_CPPOPTIONS_H */
C   ==================================================================
C-- Package-specific Options & Macros go here

C       >>> Cost function contributions
#undef ALLOW_COST_ATLANTIC
#undef ALLOW_COST_ATLANTIC_HEAT

#define ALLOW_COST_TEST
#undef ALLOW_COST_TSQUARED
#undef ALLOW_COST_TRACER

C   ==================================================================
#endif /* ndef ECCO_CPPOPTIONS_H */
#endif /* ALLOW_COST */
#endif /* COST_OPTIONS_H */
