% This script compares archetypes and cluster centers found in four
% binary datasets, both qualitatively and quantitatively
% Author: Sohan Seth, sohan.seth@hiit.fi

clear all, close all, rng default, flagSave = false;
epsilon = 10^-12; maxTrial = 1; FONT = 10;
computeLikelihood = @(x,p)(sum(x .* log(p + epsilon) + (1-x) .* log(1 - p + epsilon)));

for data = [1 2 3 4]; % 1,2 has labels, 3 does not
    
    switch data
        case 1
            load spect
            file = 'spect_save';
        case 2
            load congress
            file = 'congress_save';
        case 3
            % nSamFeat = load('bmix-1.11/data/dna_amp_chr_17.data'); [nSam, nFeat] = size(nSamFeat); label = ones(nSam, 1);
            load dna
            file = 'dna_save';
        case 4
            nSamFeat = rand(200,25) > 0.5; [nSam, nFeat] = size(nSamFeat); label = ones(nSam, 1);
            file = 'binsynth_save';
    end
    disp(file)
    
    if exist([file,'.mat']) == 0 % pre computed prototypes and associated information
        
        % -1 implies missing values, discard them from the analysis
        if any(data == [1,2])
            label(sum(features == -1,2) > 0) = [];
            features(sum(features == -1,2) > 0, :) = [];
            nSamFeat = features; [nSam, nFeat] = size(nSamFeat);
        end
        
        % Number of prototypes
        nLat = 4;
        
        % Construct input for AA
        clear nFeatSam_
        for countFeat = 1:size(nSamFeat,2)
            nFeatSam_{countFeat}(1,:) = nSamFeat(:,countFeat);
            nFeatSam_{countFeat}(2,:) = 1 - nSamFeat(:,countFeat);
        end
        
        clear matSamLat matLatSam objAA
        options = generate_options();
        for countTrial = 1:maxTrial
            [matSamLat{countTrial}, matLatSam{countTrial}, objAA{countTrial}] = paa_nominal_EM(nFeatSam_, nLat, options);
        end
        % Best AA solution indexed by indAA
        [~, indAA] = max(cellfun(@(x)(x(end)),objAA));
        % Archetypes
        matLatFeatAA = (nSamFeat' * matSamLat{indAA})';
        
        % likelihoods = [];
        % for countSam = 1:nSam
        %     for countLat = 1:nLat
        %         likelihoods(countSam, countLat) = computeLikelihood(nSamFeat(countSam, :), matLatFeatAA(countLat,:));
        %    end
        % end
        % % Cluster assignment and representative samples
        % [~, clustAA] = max(likelihoods, [], 2);
        % [~, indRepAA] = max(likelihoods);
        
        % Cluster assignment and representative samples
        [~, clustAA] = max(matLatSam{indAA});
        [~, indRepAA] = max(matLatSam{indAA}, [], 2);
        
        % Prepare data for Bernoulli mix
        clear indtmp indtmp2 obj
        dlmwrite('bmix-1.11.txt',nSamFeat,'delimiter',' ')
        for countTrial = 1:maxTrial
            [~, msg] = system(['./bmixWrapper.sh ',num2str(nFeat),' ',num2str(nLat),' bmix-1.11.txt']);
            tmp = regexp(msg,': ([+-][0-9]+.[0-9]+)','tokens');
            indtmp{countTrial} = load('small.final.clusters');
            obj{countTrial} = str2num(tmp{end}{1});
            [~, msg] = system(['tail -n ',num2str(nLat),' small.final.model']);
            indtmp2{countTrial} = str2num(msg);
        end
        % Centers and cluster assignments
        [~, ind] = max(cellfun(@(x)(x(end)),obj));
        matLatFeatBMix = indtmp2{ind};
        clustBMix = indtmp{ind};
        
        likelihoods = [];
        for countSam = 1:nSam
            for countLat = 1:nLat
                likelihoods(countSam, countLat) = computeLikelihood(nSamFeat(countSam, :), matLatFeatBMix(countLat,:));
            end
        end
        % Representative samples
        [~, indRepBMix] = max(likelihoods);
        
        %fprintf('association between AA and true %0.6f\n', nmi(clustAA,label))
        %fprintf('association between BMix and true %0.6f\n', nmi(clustBMix,label))
        %fprintf('association between AA and BMix %0.6f\n', nmi(clustBMix,clustAA))
        
        contingency = zeros(max(clustAA), max(clustBMix));
        for countSam = 1:nSam
            contingency(clustAA(countSam),clustBMix(countSam)) = ...
                contingency(clustAA(countSam),clustBMix(countSam)) + 1;
        end
        if flagSave
            save(file)
        end
        
    else
        load(file,'matLatFeatAA','matLatFeatBMix','nLat','clustAA','clustBMix','nSamFeat','nFeat','indRepAA','indRepBMix')
        % [~, indMatch] = max(contingency);
        indMatch = greedyMatch(pdist2(matLatFeatAA, matLatFeatBMix, 'cityblock'));
        if length(unique(indMatch)) ~= nLat
            fprintf('No unique matches\n')
            indMatch = [1 2 3 4];
        end
        tmp = sortrows([[1,2,3,4]', indMatch'],2); tmp = tmp(:,1);
        
        d = squareform(pdist(double(nSamFeat),'hamming'));
        fprintf('Dunns index AA: %0.6f\n',mydunns(d, clustAA))
        fprintf('Dunns index Clust: %0.6f\n',mydunns(d, clustBMix))
        fprintf('Rand index: %0.6f\n',myrands(clustAA, clustBMix))
        fprintf('Rand index: %0.6f\n',myrands(clustAA, clustBMix(randperm(length(clustBMix)))))
%       
        if 1
            
            hMain = figure('papersize',[8.5 11],'paperposition',[0 0 7 3.5]);
            colormap(flipud([[0,0,0];hot(9)]))
            
            subplot('position',[0.05 0.5 0.3 0.4]), set(gca,'fontsize',FONT),
            temp = sortrows([tmp(clustAA), nSamFeat]); ticks = find(diff(temp(:,1))); temp(:,1) = []; %temp(:,1)/max(temp(:,1));
            imagesc(temp), set(gca,'ytick',ticks,'yticklabel',{'--'}), set(gca,'xaxislocation','top') % 'xtick',[],
            set(gca,'yaxislocation','right'), text(-size(nSamFeat,2)/20, size(nSamFeat,1)/2, ['observations',char(10),'assignment to groups --'],'rotation',90,'horizontalalignment','center','interpreter','tex')
            text(1.1 * nFeat, 1, [' DI ', num2str(mydunns(d, clustAA),'%0.03f')],'fontsize',FONT,'verticalAlignment','top','edgecolor',[0.5 0.5 0.5])
            subplot('position',[0.50 0.5 0.3 0.4]), set(gca,'fontsize',FONT),
            temp = sortrows([clustBMix, nSamFeat]); ticks = find(diff(temp(:,1))); temp(:,1) = []; %temp(:,1)/max(temp(:,1));
            imagesc(temp), set(gca,'ytick', ticks, 'yticklabel',{'--'}), set(gca,'xaxislocation','top') % 'xtick',[],
            set(gca,'yaxislocation','right')
            text(1.1 * nFeat, 1, [' DI ', num2str(mydunns(d, clustBMix),'%0.03f')],'fontsize',FONT,'verticalAlignment','top','edgecolor',[0.5 0.5 0.5])
            
            h = colorbar; set(h,'position',[0.93 0.3, 0.01, 0.6],'fontsize',FONT,'xtick',[0 1]);
            
            subplot('position',[0.05 0.3 0.3 0.19]), set(gca,'fontsize',FONT)
            imagesc(matLatFeatAA(indMatch, :),[0 1]), set(gca,'ytick',[],'xtick',[]), ylabel('arch')
            subplot('position',[0.50 0.3 0.3 0.19]), set(gca,'fontsize',FONT)
            imagesc(matLatFeatBMix, [0 1]), set(gca,'ytick',[],'xtick',[]), ylabel('cen')
            
            colorLimit = max([pdist(matLatFeatAA(indMatch, :),'cityblock') / nFeat, pdist(matLatFeatBMix(indMatch, :),'cityblock') / nFeat]);
            subplot('position',[0.36 0.3 0.1 0.19]), set(gca,'fontsize',FONT)
            imagesc(squareform(pdist(matLatFeatAA(indMatch, :),'cityblock')) / nFeat / colorLimit, [0,1]), set(gca,'ytick',[],'xtick',[]), xlabel('arch')
            subplot('position',[0.81 0.3 0.1 0.19]), set(gca,'fontsize',FONT)
            imagesc(squareform(pdist(matLatFeatBMix,'cityblock')) / nFeat / colorLimit, [0,1]), set(gca,'ytick',[],'xtick',[]), xlabel('cen')
            
            subplot('position',[0.05 0.1 0.3 0.19]), set(gca,'fontsize',8)
            temp = nSamFeat(indRepAA,:);
            imagesc(temp(indMatch,:)), set(gca,'xtick',[],'ytick',[]), xlabel('features'), ylabel('arch')
            subplot('position',[0.5 0.1 0.3 0.19]), set(gca,'fontsize',8)
            temp = nSamFeat(indRepBMix,:);
            imagesc(temp), set(gca,'xtick',[],'ytick',[]), xlabel('features'), ylabel('cen')
            
            if flagSave
                saveas(gcf,[file,'.eps'],'epsc'),
                print('-dpng','-r400',[file,'.png']), 
            end
            close(hMain)
        end
    end
end