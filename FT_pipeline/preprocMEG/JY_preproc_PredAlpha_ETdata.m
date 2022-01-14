% 
% Preprocess eyetracker(ET) data for the PredAlpha project.
% 
% JY (Oct, 2019)
% 
% 


clearvars; close all; clc;

addpath(genpath('/project/3018041.02/EyeData')); %data and toolboxes
if isempty(which('ft_defaults')) %when fieldtrip is not yet in the path
    addpath('/home/common/matlab/fieldtrip');
    addpath('/home/common/matlab/fieldtrip/qsub/');
    ft_defaults;
end


SubjectIDs = 200 + [3,4,5,6,7,8,9,10,11,12,14,15,17,20,21,22,23,25,26,27,28,29,30,31,32,33,34,36,38,39,41,43,44,46];


%% read in the asc file (i.e., output of the edf2ascConverter)

% loop through subjects
for ii = 1:numel(SubjectIDs)
    
    cfg           = [];
    cfg.SubjectID = SubjectIDs(ii);
    
    % jy_batch_epochETdata( cfg );
    
    qsubfeval('jy_batch_epochETdata', cfg, ...
        'memreq',4*1024^3, 'memoverhead', 0, 'timreq', 0.8*60^2, 'timoverhead', 0, ...
        'batchid', ['ET_',num2str(cfg.SubjectID)]);
    
end