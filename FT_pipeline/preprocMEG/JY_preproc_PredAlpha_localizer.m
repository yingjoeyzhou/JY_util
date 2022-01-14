% 
% Preprocess MEG data for the PredAlpha project.
% 
% JY (July, 2019)
% 
% 


clearvars; close all; clc;

addpath(genpath('/project/3018041.02'));
if isempty(which('ft_defaults')) %when fieldtrip is not yet in the path
    addpath('/home/common/matlab/fieldtrip');
    addpath('/home/common/matlab/fieldtrip/qsub/');
    ft_defaults;
end


%% define global params for subject of interest

SubjectID = str2double( input('Specify SubjectID:\n','s') );

% define subject-specifc path and files
pathInfo = jy_definepath_predalpha( SubjectID );

% define how we wanna read our MEG data into fieldtrip
cfg0                  = [];
cfg0.dataset          = pathInfo.MEG_dsFile;
cfg0.behaviorparam    = pathInfo.BEH_paramFile;  %JY: not FT field!
cfg0.behaviordata     = pathInfo.BEH_resultsFile;%JY: not FT field!
cfg0.dosubsampling    = 1;                       %JY: not FT field!
cfg0.resamplefs       = 400;                     %JY: for ft_resampledata.
cfg0.dftfilter        = 'yes';
cfg0.dftfreq          = [50 100 150];
cfg0.padding          = 10;    
cfg0.demean           = 'no';  
cfg0.denoisesynthetic = 'no';  
cfg0.continuous       = 'yes'; 
cfg0.channel          = {'MEG'};

step1  = 'Identify and reject artefacts.';
step2  = 'Visualize the extracted components.';
answer = questdlg('Indicate the step of preprocessing:',...
    'step quest', step1, step2, step1);

switch answer
    
    case step1
        
        
        % ==== read subject's data using parameters defined by cfg0 ======
        data_loc = jy_getdata_localizerIEM( SubjectID, cfg0 );
        
        
        
        % ================== check artefact and run ICA =============
        data_loc = jy_checkartefact( pathInfo.filePreprocLog_Loc, data_loc );
        
        
        
        
        % ======== run ICA (in Donders Cluster) ============        

        ICAjob = qsubfeval('jy_batch_runICA', pathInfo.filePreprocLog_Loc, data_loc,...
                    'memreq',20*1024^3, 'memoverhead', 0, 'timreq', 8*60^2, 'timoverhead', 0, ...
                    'batchid', ['sub',num2str(SubjectID),'_ica']);
        
        
    case step2
        
        % ==== read subject's data using parameters defined by cfg0 ======
        data_loc = jy_getdata_localizerIEM( SubjectID, cfg0 );
        
        
        % ======= reject the identified artefacts for localizers ======== 
        load(pathInfo.filePreprocLog_Loc, 'toremove_overall0', 'toremove_muscle0');
        tmp        = unique([toremove_overall0(:); toremove_muscle0(:)]);
        
        cfg        = [];
        cfg.trials = setdiff( transpose(1:size(data_loc.trialinfo,1)), tmp);
        data_loc   = ft_selectdata(cfg, data_loc);
        
        
        % ======== check ICA of the data and remove a few =========
        data_loc = jy_checkcomponent( pathInfo.filePreprocLog_Loc, data_loc );
        
        
        % =========== check again the artefacts =============
        data_loc = jy_checkartefact( pathInfo.filePreprocLog_Loc, data_loc );
        
        
        % =========== save clean data ============
        save( pathInfo.filePreprocLocData, 'data_loc', '-v7.3');
        
        
        
        
        
        
        %%
        % ============= preprocess the head position data ================
        clear data_loc; %save working memory
        
        cfg0.channel = {'HLC0011','HLC0012','HLC0013', ...
            'HLC0021','HLC0022','HLC0023', ...
            'HLC0031','HLC0032','HLC0033'};
        headpos_Loc = jy_getdata_localizerIEM( SubjectID, cfg0 );
        
        % exclude trials from "checkartefact" phase 1
        load(pathInfo.filePreprocLog_Loc, 'toremove_overall0', 'toremove_muscle0');
        tmp0        = unique([toremove_overall0(:); toremove_muscle0(:)]);        
        cfg         = [];
        cfg.trials  = setdiff( transpose(1:size(headpos_Loc.trialinfo,1)), tmp0);
        headpos_Loc = ft_selectdata(cfg, headpos_Loc);
        
        % exclude trials from "checkartefact" phase 2
        load(pathInfo.filePreprocLog_Loc, 'toremove_overall', 'toremove_muscle');
        tmp         = unique([toremove_overall(:); toremove_muscle(:)]);
        cfg         = [];
        cfg.trials  = setdiff( transpose(1:size(headpos_Loc.trialinfo,1)), tmp);
        headpos_Loc = ft_selectdata(cfg, headpos_Loc);
        
        % compute displacement and rotation (rigid body transformation)
        headpos_Loc = jy_compute_headpos( headpos_Loc );
        
        % save data
        if exist(pathInfo.fileHeadpos, 'file')
            save(pathInfo.fileHeadpos, 'headpos_Loc', '-append');
        else
            save(pathInfo.fileHeadpos, 'headpos_Loc');
        end
        
        
        
        
        
        
end








