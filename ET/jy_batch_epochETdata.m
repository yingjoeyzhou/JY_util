function jy_batch_epochETdata( cfg )


pathInfo = jy_definepath_predalpha( cfg.SubjectID );


% read asc file into matlab and parse the data
if ~isempty( pathInfo.fileETasc )
    asc = read_eyelink_asc( pathInfo.fileETasc.name ); %read in the file
    [data0, event0, blinksmp, saccsmp] = asc2dat(asc); %parse asc
else
    error('Failed to find the eyelink asc file!');
end

% =========== epoch the ET data =============
% find samples of interest with fieldtrip function 
cfg                     = [];
cfg.event               = event0;
cfg.trialdef.eventtype  = 'EVENT_STIM1';
cfg.trialdef.eventvalue = 'EVENT_STIM1';
cfg.trialdef.prestim    = 1.5;
cfg.trialdef.poststim   = 3;
cfg.hdr.Fs              = asc.fsample;
cfg.hdr.nSamples        = 1 + cfg.hdr.Fs * (cfg.trialdef.poststim + cfg.trialdef.prestim);
cfg.hdr.nTrials         = numel(asc.msg);
[trl, event]            = ft_trialfun_general(cfg);

% epoch the ET data
data            = [];
data.fsample    = data0.fsample;
data.label      = data0.label;
data.sampleinfo = trl;

for ii = 1:size(trl,1)
    trlbeg = trl(ii, 1);
    trlend = trl(ii, 2);
    trlsmp = trlend - trlbeg + 1;
    data.trial{ii} = data0.trial{1}(:, trlbeg:trlend);
    data.time{ii}  = linspace( -cfg.trialdef.prestim, cfg.trialdef.poststim, trlsmp );
end

% structure it in fieldtrip friendly format
cfgTL               = [];
cfgTL.keeptrials    = 'yes';
data                = ft_timelockanalysis(cfgTL, data);
data.cfg            = [];
data.cfg.blinksmp   = blinksmp;
data.cfg.saccsmp    = saccsmp;
data.cfg.cfgPreproc = cfg;




% ======= append trialinfo to the ET data =========
% read from disk the behavioral data
load( pathInfo.BEH_resultsFile, 'mainTrl', 'primeTrl' );
load( pathInfo.BEH_paramFile, 'proc' );

% expand the trialinfo matrix: main task trials
idx = 0;
for b = 1:proc.nBlocks
    for ii = 1:(proc.nPrimeTrialsPerBlock + proc.nMainTrialsPerBlock)
        idx = idx + 1;
        
        if ii <= proc.nPrimeTrialsPerBlock
            t = ii;
            rt(idx)         = primeTrl{b}.resp(t).rt;
            correct(idx)    = primeTrl{b}.resp(t).correct;
            presence(idx)   = primeTrl{b}.type(t);
            whichTrial(idx) = 0; %JY: hard-coded
        else
            t = ii - proc.nPrimeTrialsPerBlock;
            rt(idx)         = mainTrl{b}.resp(t).rt;
            correct(idx)    = mainTrl{b}.resp(t).correct;
            presence(idx)   = mainTrl{b}.type(t);
            whichTrial(idx) = 1; %JY: hard-coded
        end
        primeType(idx) = proc.primeType( b );
        stimOri(idx)   = proc.stimOri( b );
        
        bVec(idx) = b;
        tVec(idx) = ii;
        
    end
end

% structure the info of interest
idxColumns.iBlock     = 1;
idxColumns.iTrial     = 2;
idxColumns.primeType  = 3;
idxColumns.whichTrial = 4;
idxColumns.stimOri    = 5;
idxColumns.presence   = 6;
idxColumns.rt         = 7;
idxColumns.correct    = 8;

% trialinfo matrix
trialinfo(:, idxColumns.iBlock)    = bVec;
trialinfo(:, idxColumns.iTrial)    = tVec;
trialinfo(:, idxColumns.primeType) = primeType;
trialinfo(:, idxColumns.whichTrial)= whichTrial;
trialinfo(:, idxColumns.stimOri)   = stimOri;
trialinfo(:, idxColumns.presence)  = presence;
trialinfo(:, idxColumns.rt)        = rt;
trialinfo(:, idxColumns.correct)   = correct;
data.trialinfo  = trialinfo;
data.idxColumns = idxColumns;

% ==== for each trial, check if it is excluded in the preproc. MEG data ======
data.trialinfo = horzcat( data.trialinfo, transpose( 1:size(data.trialinfo,1) ) );

% load the preprocLog of the MEG data
load(pathInfo.filePreprocLog, 'toremove_overall0', 'toremove_muscle0');
tmp0       = unique([toremove_overall0(:); toremove_muscle0(:)]);

% exclude trials from "checkartefact" phase 1
cfg        = [];
cfg.trials = setdiff( transpose(1:size(data.trialinfo,1)), tmp0);
datatmp    = ft_selectdata(cfg, data);

% exclude trials from "checkartefact" phase 2
load(pathInfo.filePreprocLog, 'toremove_overall', 'toremove_muscle');
tmp        = unique([toremove_overall(:); toremove_muscle(:)]);
cfg        = [];
cfg.trials = setdiff( transpose(1:size(datatmp.trialinfo,1)), tmp);
datatmp    = ft_selectdata(cfg, datatmp);

% mark in the trialinfo matrix which trials are excluded
excludeVec                                    = ones( size(data.trial,1), 1);
excludeVec( datatmp.trialinfo(:,end) )        = 0;
data.idxColumns.isExcluded                    = 9;
data.trialinfo(:, data.idxColumns.isExcluded) = excludeVec;


% ========== to make blink detection/rejection even easier =============
isBlink  = zeros( size(data.trial,1), numel(data.time) );
blinkBeg = data.cfg.blinksmp(:,1);
blinkEnd = data.cfg.blinksmp(:,2);
for ii = 1:size( data.trial, 1 )
    isBlink(ii,:) = data.sampleinfo(ii,1) : data.sampleinfo(ii,2);
end
for ii = 1:size(data.cfg.blinksmp,1)
    [idxR, idxC] = find( isBlink>=blinkBeg(ii) & isBlink<=blinkEnd(ii) );
    if ~isempty(idxR)
        isBlink(idxR,idxC) = 1;
    end
end
isBlink(isBlink~=1) = 0;
data.isBlink = isBlink;




% ================= save the epoched ET data ===================
save(pathInfo.filePreprocET, 'data');



end