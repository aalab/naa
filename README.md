# Archetypal Analysis for Nominal Observations

Scripts for archetypal analysis of nominal data

MATLAB scripts for reproducing figures
 - [demo](scripts/demo.m)
 - [compare archetypes and cluster centers on binary datasets](scripts/compareArchetypesClusterCenters.m)
 - [compare archetypes and number of archetypes inferred EM and VB solution respectively](scripts/compareArchetypesEMnVB.m)
 - [compare effect of hyperparameter on inferred number of archetypes](scripts/compareHyperparameterEffect.m)

C implementation
 - uses GSL CBLAS library (installation https://github.com/LuaDist/gsl), assumes path /usr/local/include/gsl/gsl_cblas.h
 - run compile.sh to generate paa_nominal_vb
 - [demo](C/demo.sh)
