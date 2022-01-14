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
cfg0.detrend          = 'no'; %JY: has to be 'no' especially for predalpha
cfg0.demean           = 'no'; %JY: has to be 'no' especially for predalpha 
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
        data_main = jy_getdata_maintask( SubjectID, cfg0 );
        
        
        
        % ================== check artefact and run ICA =============
        data_main = jy_checkartefact( pathInfo.filePreprocLog, data_main );
        
        
        
        % ======== run ICA (in Donders Cluster) ============        
        
        ICAjob = qsubfeval('jy_batch_runICA', pathInfo.filePreprocLog, data_main,...
            'memreq',60*1024^3, 'memoverhead', 0, 'timreq', 12*60^2, 'timoverhead', 0, ...
            'batchid', ['sub',num2str(SubjectID),'_ica']);
        
        
    case step2
        
        % ==== read subject's data using parameters defined by cfg0 ======
        data_main = jy_getdata_maintask( SubjectID, cfg0 );
        
        
        % ====== reject the identified artefacts for main trials =========
        load(pathInfo.filePreprocLog, 'toremove_overall0', 'toremove_muscle0');
        tmp        = unique([toremove_overall0(:); toremove_muscle0(:)]);
        
        cfg        = [];
        cfg.trials = setdiff( transpose(1:size(data_main.trialinfo,1)), tmp);
        data_main  = ft_selectdata(cfg, data_main);
        

        % ======== check ICA of the data and remove a few =========
        data_main = jy_checkcomponent( pathInfo.filePreprocLog, data_main );
        
        
        % =========== check again the artefacts =============
        data_main= jy_checkartefact( pathInfo.filePreprocLog, data_main );
        
        
        % =========== save clean data ============
        cfg         = []; 
        cfg.trials  = (data_main.trialinfo(:,end)==1);
        cfg.channel = {'MEG'};
        data        = ft_selectdata(cfg, data_main);
        data.cfg.previous = [];
        
        cfg         = [];
        cfg.channel = {'MEG'};
        cfg.trials  = (data_main.trialinfo(:,end)==0);
        cfg.latency = [-1.1, 3.8]; %JY: hard-coded
        data_prime  = ft_selectdata(cfg, data_main);
        data_prime.cfg.previous = [];
        
        save( pathInfo.filePreprocData, 'data_prime', 'data', '-v7.3'); 
        
        
        
        
        
        
        %% 
        % ============= preprocess the head position data ================
        clear data_main; %save working memory
        
        cfg0.channel = {'HLC0011','HLC0012','HLC0013', ...
            'HLC0021','HLC0022','HLC0023', ...
            'HLC0031','HLC0032','HLC0033'};
        headpos = jy_getdata_maintask( SubjectID, cfg0 );
        
        % load the preprocLog of the MEG data
        load(pathInfo.filePreprocLog, 'toremove_overall0', 'toremove_muscle0');
        tmp0       = unique([toremove_overall0(:); toremove_muscle0(:)]);
        
        % exclude trials from "checkartefact" phase 1
        cfg        = [];
        cfg.trials = setdiff( transpose(1:size(headpos.trialinfo,1)), tmp0);
        headpos    = ft_selectdata(cfg, headpos);
        
        % exclude trials from "checkartefact" phase 2
        load(pathInfo.filePreprocLog, 'toremove_overall', 'toremove_muscle');
        tmp        = unique([toremove_overall(:); toremove_muscle(:)]);
        cfg        = [];
        cfg.trials = setdiff( transpose(1:size(headpos.trialinfo,1)), tmp);
        headpos    = ft_selectdata(cfg, headpos);
        
        % compute displacement and rotation (rigid body transformation)
        headpos = jy_compute_headpos( headpos );
        
        % save data
        if exist(pathInfo.fileHeadpos, 'file')
            save(pathInfo.fileHeadpos, 'headpos', '-append');
        else
            save(pathInfo.fileHeadpos, 'headpos');
        end
        
        
        
        
        
        
end






