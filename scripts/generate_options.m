function options = generate_options()

% GENERATE_OPTIONS generates optional input parameters for PAA methods
%   [options] = generate_options(nSam) returns a structure containing
%       eps, the convergence criteria, default is 10^-6,
%       verbose, switch for textual display, default is false
%       display, switch for graphical display, default is false
%       maxIter, maximum number of iterations, default is 10000
%       robust, a value between 0 and 1, 0 for k-means, 1 for archetype, default 1
%       options.priorMatSamLat = 0.5, prior for Bayesian solution
%       options.priorMatLatSam = 0.5, prior for Bayesian solution
%       options.matFeatLat = [], archetypes, if provided then project new samples on matFeatLat
%
% copyright (c) Sohan Seth, sohan.seth@hiit.fi

options.eps = 10^-6;
options.verbose = false;
options.display = false;
options.maxIter = 10000;
options.robust = 1;
options.priorMatSamLat = 0.1;
options.priorMatLatSam = 0.3;
options.matFeatLat = [];