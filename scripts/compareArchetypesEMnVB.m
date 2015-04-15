%% Term-frequency observations
close all, clear all, rng default

MARKERSIZE = 2; LINEWIDTH = 1; FONTSIZE = 6; flagSave = false;

d = 3; K = 5; rate = 1000; n = 100; 
gammaParam = 0.5; multParam = 1000;
if d ~= 3
    fprintf('Only valid for d = 3.\n')
end
matFeatLat = convexCircle(K);
matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
matFeatSam = matFeatLat * matLatSam;
nFeatSam{1} = mnrnd(poissrnd(multParam, n, 1), matFeatSam')';

hMain = figure('papersize',[8.5 11],'paperposition',[0 0 8 2]);

options = generate_options();
options.maxIter = 1000; 
options.verbose = 0; 
options.eps = 10^-8; 
options.display = 0;

options.priorMatSamLat = 0.1; options.priorMatLatSam = 0.5;

clear objListEM objListVB VBrecord
KList = 3:10; maxTrial = 10; KPlot = 10; colorList = summer(maxTrial);
for countK = KList
    for countTrial  = 1:maxTrial
        if ~exist('selectArchetypeVB_1.mat')
            
            [matSamLatEM, matLatSamEM, obj] = paa_ordinal_EM(nFeatSam, countK, options);
            objListEM(countK == KList, countTrial) = obj(end);
            
            [matSamLatVB, matLatSamVB, obj] = paa_ordinal_VB(nFeatSam, countK, options);
            objListVB(countK == KList, countTrial) = obj(end);
            
            matSamLatEMList{countK == KList, countTrial} = matSamLatEM;
            matLatSamEMList{countK == KList, countTrial} = matLatSamEM;
            
            matSamLatVBList{countK == KList, countTrial} = matSamLatVB;
            matLatSamVBList{countK == KList, countTrial} = matLatSamVB;
            
            VBrecord(countK == KList, countTrial) = sum(max(bsxfun(@rdivide, matLatSamVB, sum(matLatSamVB)),[],2) > 0.7);
            
        else
            load('selectArchetypeVB_1.mat')
        end
        
        if countK == KPlot
            matSamLatEM = matSamLatEMList{countK == KList, countTrial};
            matLatSamEM = matLatSamEMList{countK == KList, countTrial};
            
            if countTrial == 1
                TTT = bsxfun(@rdivide, nFeatSam{1}, sum(nFeatSam{1}));
                figure(hMain)
                subplot(2,4,[1,5]),
                plot(TTT(1,:),TTT(2,:),'o','markersize',MARKERSIZE,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor',[0.5 0.5 0.5]);
                subplot(2,4,[3,7]),
                hold on
                plot(-1,-1,'marker','o','markersize',MARKERSIZE+2,'markerfacecolor',colorList(1,:),'markeredgecolor',colorList(1,:),'linestyle','none')
                plot(-1,-1,'marker','o','markersize',MARKERSIZE+2,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor',colorList(1,:),'linestyle','none')
                plot(TTT(1,:),TTT(2,:),'o','markersize',MARKERSIZE,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor',[0.5 0.5 0.5])
                HL = legend('Active','Inactive','Obs.');
                set(HL,'color','none','edgecolor',[1 1 1])
                set(gca,'fontsize',FONTSIZE)
            end
            
            subplot(2,4,[1,5]), 
            set(gca,'fontsize',FONTSIZE), box on            
            matFeatLat = TTT * matSamLatEM;
            hold on
            plot(matFeatLat(1,:),matFeatLat(2,:),'ro','markersize',MARKERSIZE+2,'markerfacecolor',colorList(countTrial,:),'markeredgecolor',colorList(countTrial,:))
            chullInd = convhull(matFeatLat(1:2,:)');
            plot(matFeatLat(1,chullInd),matFeatLat(2,chullInd),':','color',colorList(countTrial,:),'linewidth',LINEWIDTH)
            hold off, grid on, %axis([0 1 0 1 0 1])
            title('(a) Example EM solution')
            
            matSamLatVB = matSamLatVBList{countK == KList, countTrial};
            matLatSamVB = matLatSamVBList{countK == KList, countTrial};
            
            subplot(2,4,[3,7]), 
            set(gca,'fontsize',FONTSIZE), box on
            matFeatLat = TTT * bsxfun(@rdivide, matSamLatVB, sum(matSamLatVB));
            hold on
            active = max(bsxfun(@rdivide, matLatSamVB, sum(matLatSamVB)),[],2) > 0.1;
            plot(matFeatLat(1,active),matFeatLat(2,active),'ro','markersize',MARKERSIZE+2,'markerfacecolor',colorList(countTrial,:),'markeredgecolor',colorList(countTrial,:))
            plot(matFeatLat(1,~active),matFeatLat(2,~active),'ro','markersize',MARKERSIZE+2,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor',colorList(countTrial,:))
            chullInd = convhull(matFeatLat(1:2,:)');
            plot(matFeatLat(1,chullInd),matFeatLat(2,chullInd),':','color',colorList(countTrial,:),'linewidth',LINEWIDTH)
            hold off, grid on, %axis([0 1 0 1 0 1])
            title('(b) Example VB solution')
        end
        
        if countK == KPlot && countTrial == maxTrial
            subplot(2,4,6), set(gca,'fontsize',FONTSIZE)
            imagesc(matLatSamEM), colormap(flipud(gray))
            ylabel('Archetypes')
            xlabel('Samples')
            title('(a'''') Example coefficient matrix, EM')
            
            subplot(2,4,8), set(gca,'fontsize',FONTSIZE)
            imagesc(matLatSamVB), colormap(flipud(gray))
            ylabel('Archetypes')
            xlabel('Samples')
            title('(b'''') Example coefficient matrix, VB')
        end
        
        fprintf('[K = %d, Trial = %d]\n', countK, countTrial)
    end
end

figure(hMain)
subplot(2,4,2), set(gca,'fontsize',FONTSIZE), box on
hold on
plot(KList, objListEM, 'o', 'markersize', MARKERSIZE, 'markerfacecolor', [0.5 0.5 0.5],'markeredgecolor',[0.5 0.5 0.5])
plot(KList, mean(objListEM,2),'-.k','linewidth',LINEWIDTH)
ylabel('Objective function')
xlabel('Given number of archetypes K')
hold off
title('(a'') Elbow criterion')

figure(hMain)
subplot(2,4,4), set(gca,'fontsize',FONTSIZE), box on
hold on
for countK = KList
    for countK2 = KList
        if mean(VBrecord(countK == KList, :) == countK2) ~= 0
            plot(countK, countK2, 'o', 'markersize', 2 * MARKERSIZE *  mean(VBrecord(countK == KList, :) == countK2) + 1, 'markerfacecolor', [0.5 0.5 0.5],'markeredgecolor',[0.5 0.5 0.5])
        end
    end
end
ylabel('# of inferred arch.')
xlabel('Given number of archetypes K')
hold off, axis([min(KList)-1 max(KList)+1 min(KList)-1 7]) % max(KList)+1
title('(b'') Redundancy criterion')

subplot(2,4,[1,5]);
set(gca,'xlim',[0 1],'ylim',[0 0.8])
subplot(2,4,[3,7]);
set(gca,'xlim',[0 1],'ylim',[0 0.8])

if flagSave
    save selectArchetypeVB_1.mat objListEM objListVB VBrecord nFeatSam ...
        matSamLatEMList matLatSamEMList matSamLatVBList matLatSamVBList
    figureName = 'selarch.eps';
    print('-dpng','-r400',[figureName,'.png']), 
    % saveas(gcf,figureName,'epsc')
    % !epstopdf selarch.eps
end

%% Binary observations
addpath package
close all, clear all

rng default
MARKERSIZE = 2; FONTSIZE = 6;

d = 10; K = 5; n = 100; gammaParam = 0.3;
matFeatLat = rand(d,K) > 0.5; nFeat = size(matFeatLat,1);
matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
matFeatSam = rand(nFeat, n) < (matFeatLat * matLatSam);
samples = myunique(bin2dec(num2str(matFeatSam')));
archetypesTrue = sort(myunique(bin2dec(num2str(matFeatLat')),true));
%archetypesTrue = (bin2dec(num2str(matFeatLat')));

if length(archetypesTrue) ~= K
    error('Overlapping archetypes')
end

for countFeat = 1:size(matFeatLat, 1)
    nFeatSam{countFeat}(1,:) = matFeatSam(countFeat,:);
    nFeatSam{countFeat}(2,:) = 1 - nFeatSam{countFeat}(1,:);
end

hMain = figure('papersize',[8.5 11],'paperposition',[0 0 8 2]);

options = generate_options();
options.maxIter = 1000; 
options.verbose = 0; 
options.eps = 10^-6; 
options.display = 0;

options.priorMatSamLat = 0.1; options.priorMatLatSam = 0.3;

clear objListEM objListVB VBrecord
KList = 3:10; maxTrial = 10; KPlot = 8; colorList = summer(maxTrial);

for countK = KList
    for countTrial  = 1:maxTrial
        if ~exist('selectArchetypeVB_2.mat')
            [matSamLatEM, matLatSamEM, obj] = paa_ordinal_EM(nFeatSam, countK, options);
            objListEM(countK == KList, countTrial) = obj(end);
            
            [matSamLatVB, matLatSamVB, obj] = paa_ordinal_VB(nFeatSam, countK, options);
            objListVB(countK == KList, countTrial) = obj(end);
            
            VBrecord(countK == KList, countTrial) = sum(max(bsxfun(@rdivide, matLatSamVB, sum(matLatSamVB)),[],2) > 0.1);
            fprintf('[K = %d, Trial = %d]\n', countK, countTrial)
            
            matSamLatEMList{countK == KList, countTrial} = matSamLatEM;
            matLatSamEMList{countK == KList, countTrial} = matLatSamEM;
            
            matSamLatVBList{countK == KList, countTrial} = matSamLatVB;
            matLatSamVBList{countK == KList, countTrial} = matLatSamVB;
        else
            load('selectArchetypeVB_2.mat')
        end
        
        
        if countK == KPlot
            [~, maxIndEM] = max(objListEM(countK == KList,:));
            [~, maxIndVB] = max(objListVB(countK == KList,:));
            
            matSamLatEM = matSamLatEMList{countK == KList, countTrial};
            matSamLatVB = matSamLatVBList{countK == KList, countTrial};
            matLatSamEM = matLatSamEMList{countK == KList, countTrial};
            matLatSamVB = matLatSamVBList{countK == KList, countTrial};
            
            matFeatLatEM = matFeatSam * matSamLatEM > 0.5;
            matFeatLatVB = matFeatSam * bsxfun(@rdivide, matSamLatVB, sum(matSamLatVB)) > 0.5;
            %archetypesInferredEM = sort(myunique(bin2dec(num2str(matFeatLatEM')),true));
            %archetypesInferredVB = sort(myunique(bin2dec(num2str(matFeatLatVB')),true));
            archetypesInferredEM = (bin2dec(num2str(matFeatLatEM')));
            archetypesInferredVB = (bin2dec(num2str(matFeatLatVB')));
            
            activeSet = max(bsxfun(@rdivide, matLatSamVB, sum(matLatSamVB)), [], 2) > 0.1;
            
            if countTrial == 1
                archetypesTrue = sort(myunique(bin2dec(num2str(matFeatLat')),true));
                %archetypesTrue = (bin2dec(num2str(matFeatLat')));
                
                subplot(2,4,[1,5]); axis([1,countK,0,maxTrial+2]), axis off
                set(gca,'fontsize',FONTSIZE)
                title('(a) Example EM solution')  
                text(1:K,(maxTrial+1)*ones(1,K),num2str(archetypesTrue(:,1)),...
                    'fontsize',4,'backgroundcolor',[1 1 1]);%[0.6 0.6 0.6]);
                text(0,(maxTrial+1),'true',...
                    'fontsize',FONTSIZE);
                text(-0.6,5,'Trials',...
                    'fontsize',FONTSIZE,'rotation',90);
                
                subplot(2,4,[3,7]); axis([1,countK,0,maxTrial+2]), axis off
                set(gca,'fontsize',FONTSIZE)
                title('(b) Example VB solution')
                text(1:K,(maxTrial+1)*ones(1,K),num2str(archetypesTrue(:,1)),...
                    'fontsize',4,'backgroundcolor',[1 1 1]);%[0.6 0.6 0.6]);
                text(0,(maxTrial+1),'true',...
                    'fontsize',FONTSIZE);
                text(-0.5,5,'Trials',...
                    'fontsize',FONTSIZE,'rotation',90);
                
                text(1.5,-0.5,'* best','backgroundcolor',[1 1 1],'edgecolor',[1 1 1],'fontsize',FONTSIZE)
                text(3.5,-0.5,'active','backgroundcolor',[1 1 1],'edgecolor',[0 0 0],'fontsize',FONTSIZE)
                text(5.5,-0.5,'close to true','backgroundcolor',[0.8 0.8 0.8],'edgecolor',[1 1 1],'fontsize',FONTSIZE)
            end
            
            subplot(2,4,[1,5]);
            %text(1:length(archetypesInferredEM),...
                %countTrial*ones(1,length(archetypesInferredEM)),...
                %num2str(archetypesInferredEM),'fontsize',FONTSIZE,'backgroundcolor',[1 1 1]);
            [valEM, indEM] = min(pdist2(double(matFeatLatEM'),double(matFeatLat'),'Jaccard'));
            temp = zeros(1,countK); temp(indEM) = 1; indEM = logical(temp);
            text(find(indEM), maxTrial+1-countTrial*ones(1,sum(indEM)), ...
                num2str(archetypesInferredEM(indEM)),'fontsize',4,'backgroundcolor',[0.8 0.8 0.8]);
            text(find(~indEM), maxTrial+1-countTrial*ones(1,sum(~indEM)), ...
                num2str(archetypesInferredEM(~indEM)),'fontsize',4,'backgroundcolor',[1 1 1]);
            if countTrial == maxIndEM
                text(0,countTrial,[num2str(maxTrial+1-countTrial),'*'],...
                    'fontsize',FONTSIZE);
            else
                text(0,countTrial,num2str(maxTrial+1-countTrial),...
                    'fontsize',FONTSIZE);
            end
                
            subplot(2,4,[3,7]);
            %text(find(activeSet)+0.4,maxTrial+1-countTrial*ones(1,sum(activeSet))-0.01,'\_','horizontalalignment','center','edgecolor',[0 0 0])
            [valVB, indVB] = min(pdist2(double(matFeatLatVB'),double(matFeatLat'),'Jaccard'));
            temp = zeros(1,countK); temp(indVB) = 1; indVB = logical(temp);
            text(find(indVB(:) & activeSet(:)), maxTrial+1-countTrial*ones(1,sum(indVB(:) & activeSet(:))), ...
                num2str(archetypesInferredVB(indVB(:) & activeSet(:))),'fontsize',4,'backgroundcolor',[0.8 0.8 0.8],'edgecolor',[0 0 0]);
            text(find(~indVB(:) & activeSet(:)), maxTrial+1-countTrial*ones(1,sum(~indVB(:) & activeSet(:))), ...
                num2str(archetypesInferredVB(~indVB(:) & activeSet(:))),'fontsize',4,'backgroundcolor',[1 1 1],'edgecolor',[0 0 0]);
            text(find(indVB(:) & ~activeSet(:)), maxTrial+1-countTrial*ones(1,sum(indVB(:) & ~activeSet(:))), ...
                num2str(archetypesInferredVB(indVB(:) & ~activeSet(:))),'fontsize',4,'backgroundcolor',[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
            text(find(~indVB(:) & ~activeSet(:)), maxTrial+1-countTrial*ones(1,sum(~indVB(:) & ~activeSet(:))), ...
                num2str(archetypesInferredVB(~indVB(:) & ~activeSet(:))),'fontsize',4,'backgroundcolor',[1 1 1],'edgecolor',[1 1 1]);
            if countTrial == maxIndVB
                text(0,countTrial,[num2str(maxTrial+1-countTrial),'*'],...
                    'fontsize',FONTSIZE);
            else
                text(0,countTrial,num2str(maxTrial+1-countTrial),...
                    'fontsize',FONTSIZE);
            end
            
        end
        
        if countK == KPlot && countTrial == 1

            subplot(2,4,6), set(gca,'fontsize',FONTSIZE)
            imagesc(matLatSamEM), colormap(flipud(gray))
            ylabel('Archetypes')
            xlabel('Samples')
            title('(a'''') Example coefficient matrix, EM')
            
            subplot(2,4,8), set(gca,'fontsize',FONTSIZE)
            imagesc(matLatSamVB), colormap(flipud(gray))
            ylabel('Archetypes')
            xlabel('Samples')
            title('(b'''') Example coefficient matrix, VB')
        end
    end
end

figure(hMain)
subplot(2,4,2), set(gca,'fontsize',FONTSIZE), box on
hold on
plot(KList, objListEM, 'o', 'markersize', MARKERSIZE, 'markerfacecolor', [0.5 0.5 0.5],'markeredgecolor',[0.5 0.5 0.5])
plot(KList, mean(objListEM,2),'-.k')
ylabel('Objective function')
xlabel('Given number of archetypes K')
hold off
title('(a'') Elbow criterion')
set(gca,'ytick',[])

figure(hMain),
subplot(2,4,4), set(gca,'fontsize',FONTSIZE), box on
hold on
for countK = KList
    for countK2 = 1:max(KList)
        if mean(VBrecord(countK == KList, :) == countK2) ~= 0
            plot(countK, countK2, 'o', 'markersize', 2 * MARKERSIZE *  mean(VBrecord(countK == KList, :) == countK2) + 1, 'markerfacecolor', [0.5 0.5 0.5],'markeredgecolor',[0.5 0.5 0.5])
        end
    end
end
axis([0 max(KList)+1 1 max(KList)])
set(gca,'yaxislocation','right','xlim',[2 11])
%ylabel('# of inferred arch.')
text(1,0,'# of inferred arch.','fontsize',FONTSIZE,'rotation',90)
xlabel('Given number of archetypes K')
hold off, 
title('(b'') Redundancy criterion')
% 
if flagSave
    save selectArchetypeVB_2.mat objListEM objListVB VBrecord nFeatSam matFeatSam ...
        matFeatLat matSamLatEMList matLatSamEMList matSamLatVBList matLatSamVBList archetypesTrue
    figureName = 'selarch2.eps';
    print('-dpng','-r400',[figureName,'.png']), 
    % saveas(gcf,figureName,'epsc')
    % !epstopdf selarch.eps
end

%intersect(samples(:,1),archetypesInferredEM(:,1))'
%intersect(archetypesTrue(:,1),archetypesInferredEM(:,1))'
%intersect(samples(:,1),archetypesInferredVB(:,1))'
%intersect(archetypesTrue(:,1),archetypesInferredVB(:,1))'
%sort(max(bsxfun(@rdivide, matLatSamVB, sum(matLatSamVB)),[],2))

%% Bernoulli, inferring correct number of archetypes
close all, clear all, rng default
MARKERSIZE = 6; FONTSIZE = 8;

nFeat = 16; n = 1000; gammaParam = 0.3;
options = generate_options();
options.maxIter = 1000;
options.verbose = 0;
options.eps = 10^-6;
options.display = 0;
options.priorMatSamLat = 0.1;
clear VBrecord inferredK objectListVB matSamLatVBList matLatSamVBList matFeatLatList matLatSamList archetypesTrueList
paramList = 0.3; trueKList = 3:1:8; KList = 20; maxTrial = 10;

if ~exist('selectArchetypeVB_3.mat')
    for param = paramList
        options.priorMatLatSam = param;
        for countTrial = 1:maxTrial
            
            for K = trueKList
                
                archetypesTrue = 1;
                while length(unique(archetypesTrue)) ~= K
                    matFeatLat = rand(nFeat,K) > 0.5;
                    matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
                    matFeatSam = rand(nFeat, n) < (matFeatLat * matLatSam);
                    samples = myunique(bin2dec(num2str(matFeatSam')));
                    archetypesTrue = sort(myunique(bin2dec(num2str(matFeatLat')),true));
                end
                matFeatLatList{K == trueKList, param == paramList, countTrial} = matFeatLat;
                matLatSamList{K == trueKList, param == paramList, countTrial} = matLatSam;
                matFeatSamList{K == trueKList, param == paramList, countTrial} = matFeatSam;
                archetypesTrueList{K == trueKList, param == paramList, countTrial} = archetypesTrue;
                
                for countFeat = 1:size(matFeatLat, 1)
                    nFeatSam{countFeat}(1,:) = matFeatSam(countFeat,:);
                    nFeatSam{countFeat}(2,:) = 1 - nFeatSam{countFeat}(1,:);
                end
                
                for countK = KList
                    for rndTrial  = 1:10
                        
                        [matSamLatVB, matLatSamVB, obj] = paa_ordinal_VB(nFeatSam, countK, options);
                        objListVB(K == trueKList, param == paramList, countTrial, countK == KList, rndTrial) = obj(end);
                        
                        VBrecord(K == trueKList, param == paramList, countTrial, countK == KList, rndTrial) = ...
                            sum(max(bsxfun(@rdivide, matLatSamVB, sum(matLatSamVB)),[],2) > 0.1);
                        
                        matSamLatVBList{K == trueKList, param == paramList, countTrial, countK == KList, rndTrial} = matSamLatVB;
                        matLatSamVBList{K == trueKList, param == paramList, countTrial, countK == KList, rndTrial} = matLatSamVB;
                        
                        fprintf('[trueK = %d, param = %0.1f, K = %d, Trial = %d]\n', K, param, countK, countTrial)
                    end
                    
                    [~, ind] = max(squeeze(objListVB(K == trueKList, param == paramList, countTrial, countK == KList, :)));
                    inferredK(K == trueKList, param == paramList, countTrial, countK == KList) = ...
                        VBrecord(K == trueKList, param == paramList, countTrial, countK == KList, ind);
                end
            end
        end
    end
    if flagSave
        save selectArchetypeVB_3.mat objListVB VBrecord matFeatLatList matLatSamList archetypesTrueList matSamLatVBList matLatSamVBList inferredK matFeatSamList
    end
else
    load('selectArchetypeVB_3.mat')
end

hMain = figure('papersize',[8.5 11],'paperposition',[0 0 3 2]);
nList = [500, 1000];
faceColorList = [[0.5, 0.5, 0.5]; [0.3, 0.3, 0.3]];
edgeColorList = [[0.5, 0.5, 0.5]; [0.3, 0.3, 0.3]];
biasList = [-0.1 , 0.1];
hold on
plot(0, 0, 'o', 'markerfacecolor', faceColorList(1,:),'markeredgecolor',edgeColorList(1,:))
plot(0, 0, 'o', 'markerfacecolor', faceColorList(2,:),'markeredgecolor',edgeColorList(1,:))
legend({['n = ',num2str(nList(1))],['n = ',num2str(nList(2))]},'location','southeast')
for n = nList
    switch n
        case 500
            load selectArchetypeVB_3_backup.mat
        case 1000
            load selectArchetypeVB_3.mat
    end
    for param = paramList
        hold on
        temp = squeeze(inferredK(:,param == paramList,:));
        for countK = trueKList
            for countK2 = 1:max(KList)
                if mean(temp(countK == trueKList, :) == countK2) ~= 0
                    plot(countK + biasList(n == nList), countK2, 'o', 'markersize', 3 * MARKERSIZE *  ...
                        mean(temp(countK == trueKList, :) == countK2) + 1, ...
                        'markerfacecolor', faceColorList(n == nList,:),'markeredgecolor',edgeColorList(n == nList,:))
                end
            end
        end
    end
    box on, set(gca,'xlim',[2 9], 'ylim', [2 9],'fontsize',FONTSIZE)
    ylabel('# of inferred archetypes','fontsize',FONTSIZE)
    xlabel('# of true archetypes','fontsize',FONTSIZE)
end
hold off
%figureName = 'selarch3.eps';
%saveas(gcf,figureName,'epsc')

%% Simple example for paper
clear all, rng(1000),
clear J
gammaParam = 2; K = 3; n = 250;
matFeatLat{1} = [[0.6 0.3 0.1]; [0.1 0.8 0.1]; [0.2 0.0 0.8]]';
matFeatLat{2} = [[0.1 0.8 0.1]; [0.1 0.3 0.6]; [0.4 0.4 0.2]]';
matFeatLat{3} = [[0.9 0.1]; [0.3 0.7]; [0.5 0.5]]';
matFeatLat{4} = [[0.2 0.7 0.1]; [0.0 0.8 0.2]; [0.0 0.0 1]]';
matFeatLat{5} = [[0.7 0.2 0.1]; [0.0 0.2 0.8]; [0.0 0.0 1]]';
matFeatLat{6} = [[0.5 0.5]; [0.2 0.8]; [0.5 0.5]]';
matLatSam = gamrnd(gammaParam, 1, K, n); matLatSam = bsxfun(@rdivide, matLatSam, sum(matLatSam));
for count = 1:length(matFeatLat)
    matFeatSam{count} = matFeatLat{count} * matLatSam;
    %matFeatSam{count} = bsxfun(@rdivide, matFeatSam{count}, sum(matFeatSam{count}));
    nFeatSam{count} = mnrnd(1,matFeatSam{count}')';
    [J(count,:), ~] = ind2sub(size(nFeatSam{count}), find(nFeatSam{count}));
end
matLatSamTrue = matLatSam;
matFeatLatTrue = matFeatLat;
clear matLatSam

addpath package/
options = generate_options();
options.eps = 10^-6;
options.verbose = true;
options.display = false;
options.maxIter = 500;
options.priorMatSamLat = 0.1;
options.priorMatLatSam = 0.3;
EM = true;
for countRun = 1:10
    if EM
        [matSamLat{countRun}, matLatSam{countRun}, obj{countRun}] = paa_ordinal_EM(nFeatSam, 3, options);
    else
        [matSamLat{countRun}, matLatSam{countRun}, obj{countRun}] = paa_ordinal_VB(nFeatSam, 20, options);
    end
end
costs = cellfun(@(x)(x(end)),obj);
[~, ind] = max(costs);

matLatSam = matLatSam{ind};
matSamLat = matSamLat{ind};
if ~EM
    ind = max(bsxfun(@rdivide, matLatSam, sum(matLatSam)), [], 2) > 0.1;
    temp = matSamLat(:, ind);
    temp = bsxfun(@rdivide, temp, sum(temp));
else
    temp = matSamLat;
end
for count = 1:length(matFeatLat)
    matFeatLat{count} = nFeatSam{count} * temp;
end

cell2mat(matFeatLat') 
cell2mat(matFeatLatTrue')

subplot(211); imagesc(matLatSam)
subplot(212); imagesc(matLatSamTrue)

% save simpleExample.mat matFeatLat matLatSam matSamLat matLatSamTrue matFeatLatTrue nFeatSam
%bsxfun(@rdivide, matLatSam(ind,:), sum(matLatSam(ind,:)))