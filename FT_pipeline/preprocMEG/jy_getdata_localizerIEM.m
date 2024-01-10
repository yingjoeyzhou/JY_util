function data_loc = jy_getdata_localizerIEM( SubjectID, cfg0 )

% define trials
load( cfg0.behaviorparam, 'triggers', 'locProc');
load( cfg0.behaviordata, 'locTrl');

cfg_loc                     = cfg0;
cfg_loc                     = rmfield(cfg_loc, 'resamplefs');
cfg_loc.trialdef.prestim    = 0.5; 
cfg_loc.trialdef.poststim   = 1.25 + locProc.durStim; 
cfg_loc.trialdef.eventtype  = 'UPPT001';
cfg_loc.trialdef.eventvalue = triggers.safeLoc.stimOri + [1:numel(locProc.possibleOris)]; 
cfg_loc = ft_definetrial(cfg_loc);

% read data into fieldtrip structure
data_loc           = ft_preprocessing(cfg_loc);
data_loc.Subject   = SubjectID;


% add field to indicate which column of "trialinfo" corresponds to what
nDefColumns = size(data_loc.trialinfo, 2);
if nDefColumns ~= 1, warning('nDefColumns ~= 1'); end
data_loc.idxColumns.iTriggerVal = 1;
data_loc.idxColumns.iBlock      = 2;
data_loc.idxColumns.iTrial      = 3;
data_loc.idxColumns.isOddball   = 4;
data_loc.idxColumns.respCheck   = 5;
data_loc.idxColumns.orientation = 6;



% correct for PTB orientation labeling issue
ori_trigger = cfg_loc.trl(:, end);
ori_idx     = ori_trigger - triggers.safeLoc.stimOri;
ori_label   = nan(size(ori_idx));
for ii = 1:numel(ori_idx)
    ori_label(ii) = locProc.possibleOris(ori_idx(ii));
end
corrected_label = 0 - ori_label + 180;


% add trial info
data_loc.trialinfo(:, data_loc.idxColumns.iBlock) = sort(repmat( (1:locProc.nBlocks)', [locProc.nTrialsPerBlock,1])); %iBlock
data_loc.trialinfo(:, data_loc.idxColumns.iTrial) = repmat( (1:locProc.nTrialsPerBlock)', [locProc.nBlocks,1]); %iTrial
data_loc.trialinfo(:, data_loc.idxColumns.isOddball) = locProc.stimIsOdd(:); %isOddball
data_loc.trialinfo(:, data_loc.idxColumns.respCheck) = transpose([locTrl{1}.resp(:).respCheck, ...
                                    locTrl{2}.resp(:).respCheck, locTrl{3}.resp(:).respCheck]); %pressedButton
data_loc.trialinfo(:, data_loc.idxColumns.orientation) = corrected_label;



% downsample data to 400 Hz if requested
if cfg0.dosubsampling
    
    if isempty(cfg0.resamplefs) | ~isfield(cfg0,'resamplefs');
        cfg0.resamplefs = 400;
    end
    
    cfg            = [];
    cfg.resamplefs = cfg0.resamplefs;
    cfg.demean     = 'no';
    cfg.detrend    = 'no';
    data_loc = ft_resampledata(cfg, data_loc);
    % data_loc = rmfield(data_loc, 'cfg');
    
end


end