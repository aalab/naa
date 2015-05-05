% This script compares the performance of VB solution over varying hyperparameter values

flagSave = false;
if ~exist('activeArchHyper.mat')
    dList = 10; % Dimensionality of the space
    KList = 8; % Correct number of archetypes
    nList = 250; % Number of samples
    priorMatLatSamList = (0.1:0.1:1);
    maxTrial = 20;
    gammaParam = 0.4; % Concentration of Dirichlet distribution for H
    %
    options = generate_options();
    options.maxIter = 1000;
    options.verbose = 1;
    options.eps = 10^-8;
    options.display = false;
    matSamLatVB = cell(length(dList), length(KList), length(nList), 1, length(priorMatLatSamList), maxTrial);
    matLatSamVB = cell(length(dList), length(KList), length(nList), 1, length(priorMatLatSamList), maxTrial);
    objVB = cell(length(dList), length(KList), length(nList), 1, length(priorMatLatSamList), maxTrial);
    %
    for d = dList
        for n = nList
            for K = KList
                rng(countTrial)
                options.priorMatLatSam = 0.1;
                
                for countTrial = 1:maxTrial
                    matFeatLat = rand(d, K) > 0.7;
                    while any(sum(matFeatLat') == K) || any(sum(~matFeatLat') == K)
                        % avoid all 1s in row for classical solution
                        matFeatLat = rand(d, K) > 0.7;
                    end
                    matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
                    matFeatSam = rand(d, n) < (matFeatLat * matLatSam);
                    while any(sum(matFeatSam') == n)
                        % avoid all 1s in row for classical solution
                        matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
                        matFeatSam = rand(d, n) < (matFeatLat * matLatSam);
                    end
                    
                    nFeatSam = cell(d, 1);
                    for countFeat = 1:d
                        nFeatSam{countFeat}(1, :) = matFeatSam(countFeat, :);
                        nFeatSam{countFeat}(2, :) = 1 - matFeatSam(countFeat, :);
                    end
                    
                    for priorMatLatSam = priorMatLatSamList
                        rng default
                        countK = 20; % maximum number of archetypes
                        options.priorMatSamLat = priorMatLatSam;
                        [matSamLatVB{d == dList, K == KList, n == nList, 1, priorMatLatSam == priorMatLatSamList, countTrial}, ...
                            matLatSamVB{d == dList, K == KList, n == nList, 1, priorMatLatSam == priorMatLatSamList, countTrial}, ...
                            objVB{d == dList, K == KList, n == nList, 1, priorMatLatSam == priorMatLatSamList, countTrial}] = ...
                            paa_ordinal_VB(nFeatSam, countK, options);
                    end
                    fprintf('[d = %d, K = %d, n = %d, reg = %f, Trial = %d]\n', d, K, n, priorMatLatSam, countTrial)
                end
            end
        end
    end
    
    if flagSave
        save activeArchHyper
    end
    
    % % figure, hold on
    % for priorMatLatSam = priorMatLatSamList
    %     for countTrial = 1:maxTrial
    %         activeArchetypes(priorMatLatSam == priorMatLatSamList, countTrial) = ...
    %             sum(max(matLatSamVB{1, 1, 1, 1, priorMatLatSam == priorMatLatSamList, countTrial}, [], 2) > 0.15);
    %         obj(priorMatLatSam == priorMatLatSamList, countTrial) = ...
    %             cellfun(@(x)(x(end)), objVB(1, 1, 1, 1, priorMatLatSam == priorMatLatSamList, countTrial));
    %     end
    %     [~, ind] = max(obj(priorMatLatSam == priorMatLatSamList, :));
    %     % plot(activeArchetypes(priorMatLatSam == priorMatLatSamList, :), obj(priorMatLatSam == priorMatLatSamList, :), 'o')
    %     % plot(activeArchetypes(priorMatLatSam == priorMatLatSamList, ind), obj(priorMatLatSam == priorMatLatSamList, ind), 'v')
    % end

else
    
    load activeArchHyper
    figure('papersize', [8.5 11], 'paperposition', [0 0 4 3])
    hold on
    for priorMatLatSam = priorMatLatSamList
        res = myunique(activeArchetypes(priorMatLatSam == priorMatLatSamList, :));
        for count = 1:size(res, 1)
            plot(priorMatLatSam, res(count, 1), 'marker', 'o', ...
                'markerfacecolor', [0.5 0.5 0.5], ...
                'markeredgecolor', [0.5 0.5 0.5], ...
                'markersize', res(count, 2))
        end
    end
    xlabel('\alpha', 'fontsize', 12)
    ylabel('# of inferred archetypes', 'fontsize', 12)
    set(gca, 'fontsize', 12, 'xlim', [0 1.1])
    box on
    if flagSave
        saveas(gcf, 'activeArchHyper.eps', 'epsc')
        print('-dpng','-r400','activeArchHyper.png'), 
    end
end