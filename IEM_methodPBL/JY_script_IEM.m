

cfg0              = [];
cfg0.data_dir     = '3018001.01'; %data directory
cfg0.data_file    = 'preproc_localizer_data.mat'; %preprocessed data in fieldtrip structure
cfg0.orientations = 22.5:22.5:180; %orientations presented
cfg0.timeRange    = [-0.1, 1.0]; %time window to do decoder
cfg0.tau          = 1/30;        %temporal smoothing window

%% Load data
% Task: press a button when they detected a change at fixation
data_dir  = cfg0.data_dir;
data_file = cfg0.data_file;
load(fullfile(data_dir,data_file), 'dataGratingLoc');


%% Discard localiser trials with button presses (i.e., false alarms)
keepTrials                = find(dataGratingloc.trialinfo(:,6) == 0);
dataGratingloc.trialinfo  = dataGratingloc.trialinfo(keepTrials,:);
dataGratingloc.sampleinfo = dataGratingloc.sampleinfo(keepTrials,:);
dataGratingloc.trial      = dataGratingloc.trial(keepTrials,:,:); 

%dataGratingLoc.trial is a nTrial-nChannel-nTime matrix (i.e., dimord: rpt_chan_time).


%% Prepare data
numN = size(dataGratingloc.trial,1);        % Number of trials
numF = length(dataGratingloc.label);        % Number of features/channels
numT = length(dataGratingloc.time);         % Number of time points

trial = permute(dataGratingloc.trial,[2 3 1]); %dimord: chan_time_rpt
time = dataGratingloc.time;

%phi = dataGratingloc.trialinfo(:, 4) * 22.5;
phi = cfg0.orientations(dataGratingloc.trialinfo(:,4))';

clear dataGratingloc;

% Specificy time points of interest and smooth the signal
T = find(time >= cfg0.timeRange(1), 1):cfg0.timeSteps:find(time <= cfg0.timeRange(2), 1, 'last');
T = T( ...      % Only keep time points for which averaging within tau is possible
    (T > (floor(cfg0.tau/2)+1)) ...
    & (T < (numT - floor(cfg0.tau/2))) ...
);

% Update the number of time points
numT = length(T);                   

% Average and demean
X = zeros(numF, numT, numN); %dimord: chan_time_rpt
for it = 1:numT
    t = T(it);
    
    % Average within tau
    X(:, it, :) = mean(trial(:, (1:cfg0.tau) - floor(cfg0.tau/2) - 1 + t, :), 2);
    
    % Demean over trials
    mX = mean(X(:, it, :), 3);
    X(:, it, :) = X(:, it, :) - repmat(mX, [1, 1, numN]);
end

clear trial;


% JY: note line 39 - 62 can be easily done by low-pass filtering the data
% with ft_preprocessing.. but the current code does the smoothing job (not 
% in the most efficient manner though) anyways... 


%% Decoding analysis
% Prepare data into cell format and pre-compute shrinkage
X_cell = cell(numN, 1);
gamma  = zeros(numT, 1);
for in = 1:numN
    X_cell{in} = X(:, :, in); %dimord: chan_time
    
    for it = 1:numT
        [gamma(it), ~] = shrinkage_parameters(squeeze(X(:, it, :)));
    end
end
clear X;

% Train the decoder
cfg              = [];
cfg.gamma        = gamma;
cfg.numC         = 32; %JY: hard-coded, arbitrarily defined
cfg.tuning_curve = @(X) abs(cos(X)).^5;
for ii = 1:numN
    decoder{ii} = train_temporal(cfg, X_cell, phi);
end

% Apply the decoder to new data
% To reconstruct the channel-response profile, simply structure the test 
% data X as 1-by-nTrials cell structure, in which each cell is a 
% sensor-by-time matrix (skipping the part for re-structuring X_testdata).
Y = decode_temporal_generalization(X_testdata, decoder);




%% =============== sub-function ============
function [gamma, nu, S] = shrinkage_parameters(X)
% [gamma, nu] = shrinkage_parameters(X)
%
% Dimensionality of X should be Features x Repetitions

numF = size(X, 1);
numN = size(X, 2);

m = mean(X, 2);
S = cov(X');

nu = trace(S)/numF;

z = zeros(numF, numF, numN);
for n = 1:numN
    z(:, :, n) = (X(:, n) - m)*(X(:, n) - m)';
end

gamma = ...
    (numN/((numN-1)^2)) * ...
    sum(sum(var(z, [], 3))) / ...
    sum(sum((S - nu*eye(numF)).^2) ...
);

end


%% ================= sub-function ===================
function decoder = train_temporal(cfg0, X, D)

numF = size(X{1}, 1);
numT = size(X{1}, 2);
numN = length(X);

X = reshape(cell2mat(X'), [numF, numT, numN]);

decoder = cell(numT, 1); %pre-allocate mempry

traincfg.numC         = ft_getopt( cfg0, 'numC', 32 );
traincfg.tuning_curve = ft_getopt( cfg0, 'tuning_curve', @(X) abs(cos(X)).^5 );

for it = 1:numT

    traincfg.gamma = cfg0.gamma(it);
    
    % decoder{it} = feval(cfg.train_fun, cfg.train, squeeze(X(:, it, :)), D);
    decoder{it} = trialfun_pimBH( traincfg, squeeze(X(:, it, :)), D);
    
end

end


%% ================ sub-function ================
function Y = decode_temporal_generalization(X, decoder)
% 
% INPUT:
%   X      : 1-by-nTrials cell structure, each cell is a sensor-by-time matrix.
%   decoder: 1-by-nTrainingTimePoints cell structure, the weights.
% 
% OUTPUT:
%   Y: 1-by-nTrials cell structure, each cell contains nTrainingTimePoints
%       cells. e.g., Y{1}{1} corresponds to predicted channel response
%       profile of the first trial, when the weights from the first
%       training timepoint was applied to the data X.
% 
% JY

numT_train = length(decoder);

numF        = size(X{1}, 1);
numT_decode = size(X{1}, 2);
numN        = length(X);

Y = cell(numN, 1);
for in = 1:numN %loop through the trials
    
    Y{in} = cell(numT_train, 1); %preallocate memory
    
    for it_train = 1:numT_train %decoded time series based on different decoder's trained at diff. timepoints
        Y{in}{it_train} = decoder{it_train}' * X{in};
    end    
end

end
