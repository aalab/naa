function [matSamLat, matLatSam, obj] = paa_nominal_VB(nFeatSam, nLat, options, init)

% PAA_NOMINAL_VB computes archetypal patterns from ordinal observations
%   [matSamLat,matLatSam] = paa_nominal_VB(nFeatSam,nLat) returns archetypal
%   loading matrix matSamLat, and archetypal factor matrix matLatSam, from
%   nominal observations stored in cell array nFeatSam, given number
%   of archetypes nLat. Each cell is a matrix with rows as options and 
%   columns as observation. The archetypes can be computed from loading 
%   matrix, and observation matrix as matFeatSam x matSamLat. Each nominal
%   feature can be arbitrary number of options, and multiple answers for
%   each observations.
%
%   options is an optional structure specifying paramters,
%       eps, the convergence criteria, default is 10^-6, and
%       verbose, switch for textual display, default is false
%       display, switch for graphical display, default is false
%       maxIter, maximum number of iterations, default is 10000
%
%   obj is an optional vector storing the value of objective function
%
%   copyright (c) Sohan Seth, sohan.seth@hiit.fi

% rng default; fprintf('default initialization activated for debugging.\n')

if nargin < 2
    error('Observation matrix and number of archetypes must be provided');
end

if nargin == 2
    options = generate_options();
end

epsilon = 2.2204e-16;
eps = options.eps;
verbose = options.verbose;
display = options.display;
maxIter = options.maxIter;
priorMatSamLat = options.priorMatSamLat;
priorMatLatSam = options.priorMatLatSam;
saveFileName = [];
if isfield(options,'saveFileName')
    saveFileName = options.saveFileName;
    fprintf('Saving intermediate result in %s\n',options.saveFileName)
end
clear options

if display
    figureMain = figure;
    fprintf('Only two features of the first view will be displayed.\n')
end

nFeat = length(nFeatSam);
nSam = size(nFeatSam{1}, 2);

matFeatSam = cell(nFeat,1);
for countFeat = 1:nFeat
    matFeatSam{countFeat} = bsxfun(@rdivide, nFeatSam{countFeat}, sum(nFeatSam{countFeat}));
end

obj = zeros(maxIter, 1); obj(:) = Inf;

if display && nFeat > 2
    fprintf('More than two features. Only the first two will be displayed.\n');
end

% Initialization
if exist('init','var')
    matSamLat = priorMatSamLat + init.matSamLat;
    matLatSam = priorMatLatSam + init.matLatSam;
else
    % Repeated entries have same initialization
    [val, ~, tmp] = unique(cell2mat(nFeatSam)','rows');
    matSamLat = priorMatSamLat *(1 + 0.1*rand(length(val), nLat));
    matLatSam = priorMatLatSam *(1 + 0.1*rand(nLat, length(val)));
    matSamLat = matSamLat(tmp, :); matLatSam = matLatSam(:, tmp);
end


obj(1) = computeCost(nFeatSam, matSamLat, matLatSam, priorMatLatSam, priorMatSamLat);
tic
for iter = 1:maxIter
    if verbose
        if ~mod(iter, 10)
            fprintf('.')
        end
        if ~mod(iter, 100)
            fprintf(' [%d %0.6f s]\n',iter,toc);
        end
    end
    
    % Repeat before computing obj
    for rep = 1:100
        psiMatSamLat = psi(matSamLat); psiMatSamLat = exp(bsxfun(@minus, psiMatSamLat, psi(sum(matSamLat))));
        psiMatLatSam = psi(matLatSam); psiMatLatSam = exp(bsxfun(@minus, psiMatLatSam, psi(sum(matLatSam))));
        
        temp = cell(nFeat,1);
        for countFeat = 1:nFeat
            temp{countFeat} = nFeatSam{countFeat} ./ (epsilon + matFeatSam{countFeat} * psiMatSamLat * psiMatLatSam);
        end
        
        matSamLat = zeros(nSam, nLat); matLatSam = zeros(nLat, nSam);
        for countFeat = 1:nFeat
            matSamLat = matSamLat + (matFeatSam{countFeat}' * temp{countFeat} * psiMatLatSam') .* psiMatSamLat;
            matLatSam = matLatSam + (psiMatSamLat' * matFeatSam{countFeat}' * temp{countFeat}) .* psiMatLatSam;
        end
        matSamLat = matSamLat + priorMatSamLat; matLatSam = matLatSam + priorMatLatSam;
        
        if verbose && mod(rep, 10) == 0
            fprintf('.')
        end
    end

    % Convergence    
    obj(iter+1) = computeCost(nFeatSam, matSamLat, matLatSam, priorMatLatSam, priorMatSamLat);
    
    if obj(iter+1) < obj(iter)
        if verbose
            fprintf('Lower bound dropped.\n')
        end
    else
        
        if abs((obj(iter+1) - obj(iter))/obj(iter)) < eps
            if verbose
                fprintf('\nconvergence reached in %d iterations\n', iter)
            end
            break;
        end
    end
    
    if ~isempty(saveFileName)
        save(saveFileName, 'matSamLat', 'matLatSam', 'obj', 'nLat')
        if verbose; fprintf(' [Iter %d ELBO %0.6f Time %0.6fs, saving results]\n',iter,obj(iter+1),toc); end
    else
        if verbose; fprintf(' [Iter %d ELBO %0.6f Time %0.6fs]\n',iter,obj(iter+1),toc); end
    end
    
    % Display
    if display && ~mod(iter,10)
        figure(figureMain), clf, hold on,
        TTT = bsxfun(@rdivide, nFeatSam{1}, sum(nFeatSam{1}));
        plot(TTT(1,:), TTT(2,:),'o','markerfacecolor','b','markeredgecolor',[1 1 1])
        nFeatLat = TTT * bsxfun(@rdivide, matSamLat, sum(matSamLat));
        plot(nFeatLat(1,:), nFeatLat(2,:),'o','markerfacecolor','r','markeredgecolor',[1 1 1])
        hold off,
    end
end

if iter == maxIter && verbose
    fprintf('\nmaximum iteration reached\n')    
end
%fprintf('\n')
obj(isinf(obj)) = [];

function val = computeCost(nFeatSam,...
    thetaSamLat,betaLatSam,priorMatLatSam,priorMatSamLat)

% entMulti = @(p)(-sum(p(p ~= 0) .* log(p(p ~= 0))));
epsilon = 10^-16;
psiThetaSamLat = bsxfun(@minus, psi(thetaSamLat), psi(sum(thetaSamLat)));
expPsiThetaSamLat = exp(psiThetaSamLat);
psiBetaLatSam = bsxfun(@minus, psi(betaLatSam), psi(sum(betaLatSam)));
expPsiBetaLatSam = exp(psiBetaLatSam);

val = 0;
for countFeat = 1:length(nFeatSam)
    matFeatSam = bsxfun(@rdivide, nFeatSam{countFeat}, sum(nFeatSam{countFeat}));
    val = val + sum(sum( nFeatSam{countFeat} .* ...
        log(epsilon + matFeatSam * expPsiThetaSamLat * expPsiBetaLatSam)));
end

val = val + (priorMatLatSam - 1) * sum(sum(psiBetaLatSam));
val = val + (priorMatSamLat - 1) * sum(sum(psiThetaSamLat));

val = val + sum(entDir(thetaSamLat));
val = val + sum(entDir(betaLatSam));