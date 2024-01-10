function data_main = jy_getdata_maintask( SubjectID, cfg0 )


% define trials

load( cfg0.behaviorparam, 'triggers', 'stim', 'proc');
load( cfg0.behaviordata, 'mainTrl' );

cfg_main                     = cfg0;
cfg_main                     = rmfield(cfg_main, 'resamplefs');
cfg_main.trialdef.prestim    = 2; 
cfg_main.trialdef.poststim   = round( stim.dur.target + stim.dur.ISI +stim.dur.mask + stim.dur.postISI + (stim.dur.ITI(1) + stim.dur.PreFix), 1 ); 
cfg_main.trialdef.eventtype  = 'UPPT001';
cfg_main.trialdef.eventvalue = triggers.main.stimOn; 

cfg_main = ft_definetrial(cfg_main);


% check if the number of triggers equal the number of trials
nTriggers = size( cfg_main.trl, 1);
nTrials   = (proc.nPrimeTrialsPerBlock + proc.nMainTrialsPerBlock) * proc.nBlocks;
if nTriggers == nTrials
    disp('The number of trial triggers recorded equal the number of trials programmed!');
else
    error('The number of trial triggers recorded does NOT equal the number of trials programmed!');
end

% trialtype vector: maintrials = 1; primingtrials = 0.
trialtype = ones(nTrials, 1);
beg_prime = 1 + (proc.nPrimeTrialsPerBlock + proc.nMainTrialsPerBlock)*[0:(proc.nBlocks-1)];
end_prime = proc.nPrimeTrialsPerBlock + (proc.nPrimeTrialsPerBlock + proc.nMainTrialsPerBlock)*[0:(proc.nBlocks-1)];
for b = 1:proc.nBlocks
    trialtype( beg_prime(b):end_prime(b) ) = 0;
end

% index vector
idx = 0;
for b = 1:proc.nBlocks
    for ii = 1:(proc.nPrimeTrialsPerBlock + proc.nMainTrialsPerBlock)
        idx = idx + 1;
        trialIndex(idx) = ii;
        blockIndex(idx) = b;
    end
end

% read data into fieldtrip structure
if isfield(cfg0,'channel')
    cfg_main.channel = cfg0.channel;
end
data_main           = ft_preprocessing(cfg_main);
data_main.Subject   = SubjectID;


% structure my trialinfo matrix better
data_main.idxColumns.iTriggerVal = 1;
data_main.idxColumns.iBlock      = 2;
data_main.idxColumns.iTrial      = 3;
data_main.idxColumns.iTrialType  = 4;

data_main.trialinfo(:, data_main.idxColumns.iBlock)     = blockIndex(:);
data_main.trialinfo(:, data_main.idxColumns.iTrial)     = trialIndex(:);
data_main.trialinfo(:, data_main.idxColumns.iTrialType) = trialtype(:);



% downsample data to 400 Hz if requested
if cfg0.dosubsampling
    
    if isempty(cfg0.resamplefs) | ~isfield(cfg0,'resamplefs');
        cfg0.resamplefs = 400;
    end
    
    cfg            = [];
    cfg.resamplefs = cfg0.resamplefs;
    cfg.demean     = 'no';
    cfg.detrend    = 'no';
    data_main = ft_resampledata(cfg, data_main);
    data_main = rmfield(data_main, 'cfg');
    
end



end