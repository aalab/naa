function [matSamLat, matLatSam, obj] = paa_nominal_EM(nFeatSam, nLat, options)

% PAA_NOMINAL_EM computes archetypal patterns from ordinal observations
%   [matSamLat,matLatSam] = paa_nominal_EM(nFeatSam,nLat) returns archetypal
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
%       matFeatLat, if not empty then matSamLat is not learned, but 
%               archetypes are fixed to options.matSamLat
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

eps = options.eps;
verbose = options.verbose;
display = options.display;
maxIter = options.maxIter;
matFeatLat = options.matFeatLat;
nFeatSam = nFeatSam(:);
clear options

if display
    figureMain = figure;
    fprintf('Only two features of the first view will be displayed.\n')
end

obj = zeros(maxIter + 1, 1); obj(:) = Inf;

nFeat = length(nFeatSam);
nSam = size(nFeatSam{1}, 2);

matFeatSam = cell(nFeat,1);
for countFeat = 1:nFeat
    matFeatSam{countFeat} = bsxfun(@rdivide, nFeatSam{countFeat}, sum(nFeatSam{countFeat}));
end

% Repeated entries have same initialization
[val, ~, tmp] = unique(cell2mat(nFeatSam)','rows');

% Initialization
if isempty(matFeatLat)
    matSamLat = rand(length(val), nLat); matSamLat = matSamLat(tmp, :);
    matSamLat = bsxfun(@rdivide, matSamLat, sum(matSamLat));
else
    nLat = size(matFeatLat{1}, 2);
end
matLatSam = rand(nLat, length(val)); matLatSam = matLatSam(:, tmp);
matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));

if isempty(matFeatLat)
    computeCost = @(nFeatSam, matFeatSam, matSamLat, matLatSam, epsilon)...
        (sum(sum(nFeatSam .* log(epsilon + matFeatSam * matSamLat * matLatSam),2)));
else
    computeCost = @(nFeatSam, matFeatLat, matLatSam, epsilon)...
        (sum(sum(nFeatSam .* log(epsilon + matFeatLat * matLatSam),2)));
end

epsilon = 2.2204e-16;  % Avoiding log(0) and division by 0

obj(1) = 0;
for countFeat = 1:nFeat
    if isempty(matFeatLat)
        obj(1) = obj(1) + computeCost(nFeatSam{countFeat}, matFeatSam{countFeat}, matSamLat, matLatSam, epsilon);
    else
        obj(1) = obj(1) + computeCost(nFeatSam{countFeat}, matFeatLat{countFeat}, matLatSam, epsilon);
    end
end

temp = cell(nFeat,1);
for iter = 1:maxIter
    if verbose
        if ~mod(iter, 10)
            fprintf('.')
        end
        if ~mod(iter, 100)
            fprintf(' [%d]\n',iter);
        end
    end
    
    % Repeat before computing obj
    for rep = 1:100
        % Expectation
        for countFeat = 1:nFeat
            if isempty(matFeatLat)
                temp{countFeat} = nFeatSam{countFeat} ./ ((epsilon + matFeatSam{countFeat} * matSamLat * matLatSam));
            else
                temp{countFeat} = nFeatSam{countFeat} ./ ((epsilon + matFeatLat{countFeat} * matLatSam));
            end
        end
        
        if isempty(matFeatLat)
            matSamLatNew = zeros(nSam, nLat); 
        end
        matLatSamNew = zeros(nLat, nSam);
        for countFeat = 1:nFeat
            if isempty(matFeatLat)
                matSamLatNew = matSamLatNew + (matFeatSam{countFeat}' * temp{countFeat} * matLatSam') .* matSamLat;
            end
            if isempty(matFeatLat)
                matLatSamNew = matLatSamNew + (matSamLat' * matFeatSam{countFeat}' * temp{countFeat}) .* matLatSam;
            else
                matLatSamNew = matLatSamNew + (matFeatLat{countFeat}' * temp{countFeat}) .* matLatSam;
            end
        end
        if isempty(matFeatLat)
            matSamLatNew = matSamLatNew + epsilon; 
        end
        matLatSamNew = matLatSamNew + epsilon;
        
        % Maximization
        if isempty(matFeatLat)
            matSamLat = bsxfun(@rdivide,matSamLatNew,sum(matSamLatNew));
        end
        matLatSam = bsxfun(@rdivide,matLatSamNew,sum(matLatSamNew));
        
        if verbose && mod(rep, 10) == 0
            fprintf('.')
        end
    end
    if verbose; fprintf('%d\n',iter); end
    
    % Convergence
    obj(iter+1) = 0;
    for countFeat = 1:nFeat
        if isempty(matFeatLat)
            obj(iter+1) = obj(iter+1) + computeCost(nFeatSam{countFeat}, matFeatSam{countFeat}, matSamLat, matLatSam, epsilon);    
        else
            obj(iter+1) = obj(iter+1) + computeCost(nFeatSam{countFeat}, matFeatLat{countFeat}, matLatSam, epsilon);    
        end
    end
    if abs((obj(iter+1) - obj(iter))/obj(iter)) < eps
        if verbose
            fprintf('\nconvergence reached in %d iterations\n', iter)
        end
        break;
    end
    
    % Display
    if display && ~mod(iter,10)
        figure(figureMain), clf, hold on,
        TTT = bsxfun(@rdivide, nFeatSam{1}, sum(nFeatSam{1}));
        plot(TTT(1,:), TTT(2,:),'o','markerfacecolor','b','markeredgecolor',[1 1 1])
        nFeatLat = TTT * matSamLat;
        plot(nFeatLat(1,:), nFeatLat(2,:),'o','markerfacecolor','r','markeredgecolor',[1 1 1])
        hold off,
    end
end

if iter == maxIter && verbose
    fprintf('\nmaximum iteration reached\n')    
end

obj(isinf(obj)) = [];

if ~isempty(matFeatLat)
    matSamLat = [];
end