%% Demo: archetypal analysis for nominal observations

clear all, close all
rng default,
options = generate_options(); 
options.verbose = true; 
options.display = false;
options.eps = 10^-6;
options.maxIter = 1000;
options.matFeatLat = [];
n = 500;

% Generating data
nFeat = 2; % number of features (views)
nCat = floor(rand(nFeat, 1) * 3) + 2; % number of categories in each view (2-4 for simplicity)
matFeatSam = ceil(rand(nFeat, n) .* (nCat * ones(1, n))); % for each view each entry 1-nCat

% Preparing data for paa, one cell entry for each view
nFeatSam = cell(nFeat, 1);
for count = 1:nFeat
    % each column has one 1 and rest are 0s
    nFeatSam{count} = bsxfun(@eq, matFeatSam(count, :), (1:nCat(count))');
end

% Learning archetypes using VB, maximum number of archetypes
options.priorMatLatSam = 0.3;
options.priorMatSamLat = 0.1;
[matSamLat_VB, matLatSam_VB, ~] = paa_nominal_VB(nFeatSam, 20, options);
matSamLat_VB = bsxfun(@rdivide, matSamLat_VB, sum(matSamLat_VB));
matLatSam_VB = bsxfun(@rdivide, matLatSam_VB, sum(matLatSam_VB));

% Finding active archetypes
activeArchetypes = find(max(matLatSam_VB, [], 2) > 0.15);
fprintf('number of active archetypes %d\n', length(activeArchetypes))
matSamLat_VB = matSamLat_VB(:, activeArchetypes);
matLatSam_VB = matLatSam_VB(activeArchetypes, :);
matLatSam_VB = bsxfun(@rdivide, matLatSam_VB, sum(matLatSam_VB));

% Archetypes for each view
archetypes_VB = cell(nFeat, 1);
for count = 1:nFeat
    archetypes_VB{count} = nFeatSam{count} * matSamLat_VB;
end

% Learning archetypes using EM
[matSamLat_EM, matLatSam_EM, ~] = paa_nominal_EM(nFeatSam, ...
    length(activeArchetypes), options);
% Archetypes for each view
archetypes_EM = cell(nFeat, 1);
for count = 1:nFeat
    archetypes_EM{count} = nFeatSam{count} * matSamLat_EM;
end

% Computing projections given archetypes
options.matFeatLat = archetypes_EM;
options.display = false;
[~, matLatSam, ~] = paa_nominal_EM(nFeatSam, [], options);
options.matFeatLat = [];

for count = 1:nFeat
fprintf('difference between projections %0.6f\n', ...
    norm(archetypes_EM{count} * matLatSam - archetypes_EM{count} * matLatSam_EM, 'inf') ...
    / numel(nFeatSam{count}))
end