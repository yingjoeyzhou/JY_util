%
% Figure N: AUC of the loc2task and prime2task (overall) generalization.
% 
%  * loc2task AUC tem. gen. matrix.
%  * prime2task AUC tem. gen. matrix.
%  * diagonal of the loc2task performance: conservative vs. liberal.
%  * diagonal of the prime2task performance: conservative vs. liberal.
%
% JY (Apr, 2020)
% 

clearvars; close all; clc;

addpath(genpath('/project/3018041.02'));
addpath(genpath('/home/predatt/yinzho/Documents/cbrewer'));
addpath(genpath('/home/predatt/yinzho/Documents/shadedErrorBar'));
addpath(genpath('/home/predatt/yinzho/Documents/export_fig'));
if isempty(which('ft_defaults')) %when fieldtrip is not yet in the path
    addpath('/home/common/matlab/fieldtrip');
    addpath('/home/common/matlab/fieldtrip/qsub/');
    ft_defaults;
end

% subject list
SubjectIDs = readSubjectIDs_manuscript;

% output .mat file
outputmat = 'Figure_decoding0.mat';

% conditions of interests
PRMval  = [2, 8];
PRM     = {'Conservative','Liberal'};
PRESval = [1, 0];
PRES    = {'Present','Absent'};

% set figure defaults
printPNG = 1;



%% loop through subjects to fetch data
for iSub = 1:numel(SubjectIDs)
    
    
    fprintf('\n\n starting to process subject #%d...\n',iSub);
    
    % ===================== prime2task ====================
    f_p2m = sprintf('prime2task_PA_logreg_%s.mat', num2str(SubjectIDs(iSub)) );
    load(f_p2m, 'cfg0', 'classifier_dval', 'classifier_prob', 'classifier_resp', 'perf_mat','auc_mat');
    presence = cfg0.label_test; 
    prime    = cfg0.trialinfo_test(:, cfg0.idxColumns_test.primeType);
    corr     = cfg0.trialinfo_test(:, cfg0.idxColumns_test.correct); %subject's correctness
    correct  = ( classifier_resp - permute( presence, [3,2,1] ) == 0 ); %decoder's correctness
    
    % overall accuracy
    auc_p2m(:,:,iSub) = auc_mat; %computed but not shown
    acc_p2m(:,:,iSub) = perf_mat;
    p_p2m.cfg0        = cfg0;
    
    
    % accuracy by condition
    for iPrm = 1:2
        sel = ( prime == PRMval(iPrm) );
        acc_p2m_Prm{iPrm}(:,:,iSub) = mean( correct(:,:,sel), 3 );
    end
    
    
    % accuracy by condition: Hit vs. Miss
    acc_p2m_HM{1}(:,:,iSub) = mean( correct(:,:,(presence==2 & corr~=0)), 3 ); %hit trials
    acc_p2m_HM{2}(:,:,iSub) = mean( correct(:,:,(presence==2 & corr==0)), 3 ); %miss trials
    
    acc_p2m_CF{1}(:,:,iSub) = mean( correct(:,:,(presence==1 & corr~=0)), 3 ); %CR trials
    acc_p2m_CF{2}(:,:,iSub) = mean( correct(:,:,(presence==1 & corr==0)), 3 ); %FA trials
    
    
    % ====================== loc2task ============================
    f        = sprintf('loc2task_CW_CCW_logreg_%s.mat',num2str(SubjectIDs(iSub)));
    out      = load( f, 'cfg0', 'resp_decoder','dval_decoder');
    p        = jy_definepath_predalpha( out.cfg0.SubjectID );
    [ti, iC] = jy_definetrialinfo_predalpha( out.cfg0.trialinfo,  out.cfg0.idxColumns, p, 'main' );
    presence = ti(:, iC.presence );
    ori      = out.cfg0.label_test; %ti(:, iC.ori );
    prime    = ti(:, iC.primeType );
    
    
    % overall accuracy
    correct = (out.resp_decoder - permute( ori, [3,2,1] ) == 0);
    acc_l2m(:,:,iSub) = mean( correct, 3);
    acc_l2m_pres(:,:,iSub) = mean( correct(:,:,presence==1),3 );
    acc_l2m_abs(:,:,iSub)  = mean( correct(:,:,presence==0),3 );
    
    % accuracy by condition: presence x priming type
    for iPres = 1:2
        for iPrm = 1:2
            sel = ( presence == PRESval(iPres) ) & ( prime == PRMval(iPrm) );
            acc_l2m_PresxPrm{iPres,iPrm}(:,:,iSub) = mean( correct(:,:,sel), 3 );
        end
    end
    
    %{
        for iRow = 1:size(out.dval_decoder,1)
            for iCol = 1:size(out.dval_decoder,2)
                
                for iPres = 1:2
                    sel = ( presence == PRESval(iPres) );
                    auc_l2m{iPres}(iRow,iCol,iSub) = mv_calculate_performance('auc', 'dval', squeeze(out.dval_decoder(iRow,iCol,sel)), ori(sel));
                end
                
                if iRow==iCol %compute the diagonal
                    for iPres = 1:2
                        for iPrm = 1:2
                            sel = ( prime == PRMval(iPrm) ) & (presence == PRESval(iPres));
                            diag_auc_l2m{iPres}(iRow, iPrm, iSub) = mv_calculate_performance('auc', 'dval', squeeze(out.dval_decoder(iRow,iCol,sel)), ori(sel));
                        end
                    end
                end
                
            end
        end
    %}
    p_l2m.cfg0 = out.cfg0;
    
    
end




%% Do stats


% find the peak of decoding
diagACC = diag( mean( acc_p2m , 3) );
fprintf( '\n\n\nPrime2Main decoding reaches the peak at time = %s sec.\n', num2str(time_p2m( diagACC == max( diagACC ) )) );


% time of the full matrix
time_p2m = linspace( p_p2m.cfg0.latency(1), p_p2m.cfg0.latency(2), size(acc_p2m,1));
time_l2m = linspace( p_l2m.cfg0.latency(1), p_l2m.cfg0.latency(2), size(acc_l2m,1));

if exist( outputmat, 'file')
    
    load( outputmat, 's_p2m', 'pm_p2m', 'nm_p2m', ...
                's_l2m', 'pm_l2m', 'nm_l2m',...
                's_l2m_P', 'pm_l2m_P', 'nm_l2m_P');
    
else %run stat
    
    % compare to chance-level decoding (i.e., ACC = 0.5)
    cfgIN                  = [];
    cfgIN.clustertail      = 0;
    cfgIN.clusteralpha     = 0.05;
    cfgIN.tail             = 0;
    cfgIN.alpha            = 0.025;
    cfgIN.numrandomization = 10000;
    cfgIN.clusterstatistic = 'maxsum';
    
    cfgIN.traintime         = time_p2m;
    cfgIN.testtime          = time_p2m;
    cfgIN.stattimewin       = [0, 1.0];
    [s_p2m, pm_p2m, nm_p2m] = jy_permutationtest4tempgenmat(cfgIN, acc_p2m, 0.5*ones( size(acc_p2m) ) ); %compare to chance
    
    cfgIN.traintime         = time_l2m;
    cfgIN.testtime          = time_l2m;
    cfgIN.stattimewin       = [0, 0.3];
    [s_l2m, pm_l2m, nm_l2m] = jy_permutationtest4tempgenmat(cfgIN, acc_l2m, 0.5*ones( size(acc_l2m) ) ); %compare to chance
    
    [s_l2m_P, pm_l2m_P, nm_l2m_P] = jy_permutationtest4tempgenmat(cfgIN, acc_l2m_pres, 0.5*ones( size(acc_l2m_pres) ) ); %compare to chance
    
    
    save( outputmat, 's_p2m', 'pm_p2m', 'nm_p2m', ...
        's_l2m', 'pm_l2m', 'nm_l2m',...
        's_l2m_P', 'pm_l2m_P', 'nm_l2m_P');
end



%% temgen Prime2Main

figcfg        = [];
figcfg.zlim   = 0.5 + [-0.2, 0.2]; 
figcfg.xlim   = [-0.2, 1.0];
figcfg.xticks = [0:0.5:1.0];
figcfg.ylim   = [-0.2, 1.0];
figcfg.yticks = [0:0.5:1.0];
figcfg.CMAP   = flipud( colorcet( 'd13' ) );

figure(1), clf, hold on,
set( gcf, 'units', 'centimeters', 'position', [0 0 6 6]);
plotDecodingTimeGenMat( figcfg, time_p2m, acc_p2m, pm_p2m+nm_p2m );
% JY: export_fig cannot handle surfaceplots with contours on it and export
% .svg format as expected. hence, we asked for .png as output format.
if printPNG, export_fig( 'decoding0_p2m', '-PNG', '-r500','-nocrop', '-transparent'); end



%% colorbar
% JY: stupidly, we need to export the colorbar separately
figure(1), clf, hold on,
imagesc( randn(300,400), figcfg.zlim ); colormap( figcfg.CMAP ); colorbar;
if printPNG, print( 'colorbar_decoding0_p2m', '-dsvg'); end




%% temgen Loc2Main

figcfg        = [];
figcfg.zlim   = 0.5 + [-0.1, 0.1]; 
figcfg.xlim   = [-0.1, 1.0];
figcfg.xticks = [0:0.5:1.0];
figcfg.ylim   = [-0.1, 1.0];
figcfg.yticks = [0:0.5:1.0];
figcfg.CMAP   = flipud( colorcet( 'd13' ) );

figure(201), clf, hold on,
set( gcf, 'units', 'centimeters', 'position', [0 0 9 6]);
plotDecodingTimeGenMat( figcfg, time_l2m, acc_l2m, pm_l2m );
if printPNG, export_fig( 'decoding0_l2m_alltrials', '-PNG', '-r500','-nocrop', '-transparent'); end


figure(202), clf, hold on,
set( gcf, 'units', 'centimeters', 'position', [0 0 9 6]);
plotDecodingTimeGenMat( figcfg, time_l2m, squeeze( mean( cat(4, acc_l2m_PresxPrm{1,:}), 4)), pm_l2m_P );
if printPNG, export_fig( 'decoding0_l2m_presenttrials', '-PNG', '-r500','-nocrop', '-transparent'); end



%% diagonals: prime2main
mat_acc_p2m = cat(4, acc_p2m_Prm{:});
mat_acc_l2m = cat(4, acc_l2m_PresxPrm{1,:});
for iSub = 1:size(acc_p2m,3)
    for iPrm = 1:2
        diag_acc_p2m(:,iPrm,iSub) = diag( mat_acc_p2m(:,:,iSub,iPrm) );
        diag_acc_l2m(:,iPrm,iSub) = diag( mat_acc_l2m(:,:,iSub,iPrm) );
    end
end

figure(301), clf, hold on,
set( gcf, 'units', 'centimeters', 'position', [0 0 8 6]);
cmap = cbrewer( 'div','PRGn',5, 'pchip');
cmap = cmap( [1,end],:);
figcfg        = [];
figcfg.xlim   = [-0.2, 1.0];
figcfg.xticks = [-0.2:0.2:1.0];
figcfg.ylim   = [0.45, 0.7];
figcfg.yticks = [0.45:0.05:0.7];
figcfg.cmap   = cmap;

plotDecodingTimeSeries( figcfg, time_p2m, diag_acc_p2m, []  );
if printPNG, export_fig( 'decoding0_p2m_diag_alltrials', '-svg', '-r500', '-nocrop', '-transparent'); end





%% Diagonal of Hit vs. Miss; CR vs. FA
mat_acc_HM  = cat(4, acc_p2m_HM{:});
mat_acc_CF  = cat(4, acc_p2m_CF{:});
for iSub = 1:size(acc_p2m,3)
    for iCorr = 1:2 %target type: present vs. absent
        diag_acc_p2m_HM(:,iCorr,iSub) = diag( mat_acc_HM(:,:,iSub,iCorr) );
        diag_acc_p2m_CF(:,iCorr,iSub) = diag( mat_acc_CF(:,:,iSub,iCorr) );
    end
end

prev = load( 'Figure_alphasourcetmap.mat', 'subjCF','subjHM','SubjectIDs');
SubjectCF = prev.SubjectIDs( prev.subjCF );
SubjectHM = prev.SubjectIDs( prev.subjHM );

figcfg        = [];
figcfg.xlim   = [-0.2, 1.0];
figcfg.xticks = [-0.2:0.2:1.0];
figcfg.ylim   = [0.45, 0.7];
figcfg.yticks = [0.45:0.05:0.7];
figcfg.cmap   = [0 0 0; 0.7 0.7 0.7]; %black and gray

figure(401), clf, hold on,
set( gcf, 'units', 'centimeters', 'position', [0 0 8 6]);
[~, selSubject, ~] = intersect(SubjectIDs, SubjectHM);
plotDecodingTimeSeries( figcfg, time_p2m, diag_acc_p2m_HM(:,:,selSubject), (p1 | n1)  );
%if printPNG, export_fig( 'decoding0_p2m_diag_HM', '-svg', '-r500', '-nocrop', '-transparent'); end

figure(501), clf, hold on,
iPeak = find( abs(time_p2m - 0.12) < 1e-8 ); %JY: hard-coded 0.12
bar( 1, mean( diag_acc_p2m_HM(iPeak,1,selSubject), 3), 'edgecolor',figcfg.cmap(1,:), 'facecolor',figcfg.cmap(1,:) );
bar( 2, mean( diag_acc_p2m_HM(iPeak,2,selSubject), 3), 'edgecolor',figcfg.cmap(2,:), 'facecolor',figcfg.cmap(2,:) );
ylim( figcfg.ylim );
xlim( [-0.5, 0.5] + [1, 2] );
xticks( [1,2] );

matHM = [ permute( diag_acc_p2m_HM(iPeak,1,:), [3,1,2]), permute( diag_acc_p2m_HM(iPeak,2,:), [3,1,2])];
matCF = [ permute( diag_acc_p2m_CF(iPeak,1,:), [3,1,2]), permute( diag_acc_p2m_CF(iPeak,2,:), [3,1,2])];
simple_mixed_anova( cat(3, matHM,matCF), [], {'correctness','presence'} )


%%
figure(402), clf, hold on,
set( gcf, 'units', 'centimeters', 'position', [0 0 8 6]);
[~, selSubject, ~] = intersect(SubjectIDs, SubjectCF);
plotDecodingTimeSeries( figcfg, time_p2m, diag_acc_p2m_CF(:,:,selSubject), ( p2 | n2)  );
%if printPNG, export_fig( 'decoding0_p2m_diag_CF', '-svg', '-r500', '-nocrop', '-transparent'); end









%% ==================== sub-function =====================
function plotDecodingTimeGenMat( figcfg, timevec, data, mask )
% data is a traintime-by-testtime(-by-subject) matrix.
% mask is a traintime-by-testtime logical matrix.

zmax = figcfg.zlim(2);

if size(data,3)>1
    data = mean(data, 3); 
end


uimagesc( timevec, timevec, data, figcfg.zlim); hold on,
colormap(figcfg.CMAP); 
% colorbar;
axis xy; 
axis square;

plot3( zeros(1,2), [timevec(1),timevec(end)], ones(1,2)*zmax, 'k:', 'linewidth',1);
plot3( [timevec(1),timevec(end)], zeros(1,2), ones(1,2)*zmax, 'k:', 'linewidth',1);
plot3( [timevec(1),timevec(end)], [timevec(1),timevec(end)], ones(1,2)*zmax, 'k:', 'linewidth',1); %diag

if sum( mask,'all') > 1
ct           = jy_plot_contourmask( timevec, timevec, mask );
ct.LineWidth = 0.5; %1
ct.Color     = 'k';
end

% zlim( figcfg.zlim ); 
ylim( figcfg.ylim ); 
yticks( figcfg.yticks );
ylabel( 'training time' );
xlim( figcfg.xlim ); 
xticks( figcfg.xticks );
xlabel('testing time');
set(gca,'FontName','Arial','FontSize',10,'Layer','top','TickDir','out','YColor','k','XColor','k','LineWidth',0.75);


end



%% =================== sub-function ===================
function plotDecodingTimeSeries( figcfg, timevec, data, maskvec )
% data is a time-by-condition-by-subject matrix.

mData  = squeeze( nanmean( data, 3 ) );
sdData = std( data - mean( data, [1,2]) + mean(data(:)), [], 3 ) ./ sqrt( size(data,3) );


for ii = 1:size(data,2) %loop through conditions
    shadedErrorBar( timevec, mData(:,ii), sdData(:,ii), 'lineprop',{'color',figcfg.cmap(ii,:), 'linewidth',2}, 'transparent',1  );
    hold on,
end

% denote significance
plot( timevec(logical(maskvec)), repmat( figcfg.ylim(1) + 0.1*range(figcfg.ylim), [1,sum(maskvec)] ), ...
    '-', 'color',repmat(0.6,[1,3]), 'linewidth',1.5 );

ylim( figcfg.ylim );
yticks( figcfg.yticks );
xlim( figcfg.xlim ); 
xticks( figcfg.xticks );
xlabel('Time (s)');
set(gca,'FontName','Arial','Layer','top','TickDir','out','YColor','k','XColor','k','LineWidth',1);


end


