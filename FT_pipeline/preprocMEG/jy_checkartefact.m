function data = jy_checkartefact( art_filename, data )

if ~exist(art_filename) %no such a file at all

    currStep = 'before ICA';
    
    var1 = 'toremove_overall0';
    var2 = 'toremove_muscle0';
    
else %there's a file with identified artefacts
    
    currStep = 'after ICA';
    
    var1 = 'toremove_overall';
    var2 = 'toremove_muscle';
    
end

fprintf('\n\n Checking artefacts (%s)!\n\n',currStep); 


data.trialinfo = [data.trialinfo, transpose(1:size(data.trialinfo,1)) ];


%% visual-reject overall artefacts
% Note: only delete the trials/channels that're unreasonably crazy.
% visual reject crazy trials
cfg                   = [];
cfg.renderer          = 'painters';
cfg.channel           = 'MEG';
tmp_data_overall      = ft_rejectvisual(cfg, data);
[~, toremove_overall0]= setdiff( data.trialinfo, tmp_data_overall.trialinfo, 'row');

clear tmp_data_overall
fprintf('\n\n ft_rejectvisual finished on examining overall artifacts!\n\n');


%% visual-reject muscle artefacts 
% Note: only delete the trials/channels that're unreasonably crazy.
% fill the 'nan' with zero so that fieldtrip can apply its filter
iNan = cell(size(data.trial));
for ii = 1:length(data.trial)
    iNan{ii} = isnan(data.trial{ii});
    data.trial{ii}(iNan{ii}) = 0;
end

% apply highpass filter
cfg               = [];
cfg.hpfilter      = 'yes';
cfg.hpfreq        = 100;
cfg.channel       = 'MEG';
tmp_data_filtered = ft_preprocessing(cfg, data);

% fill the 'nan's back to where they were supposed to be
for ii = 1:numel(data.trial)
    tmp_data_filtered.trial{ii}(iNan{ii}) = nan;
end

% run fieldtrip rejectvisual function
cfg                 = [];
cfg.renderer        = 'painters';
cfg.channel         = 'MEG';
tmp_data_muscle     = ft_rejectvisual(cfg, tmp_data_filtered);
[~,toremove_muscle0]= setdiff( data.trialinfo, tmp_data_muscle.trialinfo, 'row');
clear tmp_data_muscle
    
    

%% remove these crazy trials
tmp = unique([toremove_overall0(:); toremove_muscle0(:)]);
fprintf('\n\n %d trials  in total will be rejected. \n\n', numel(tmp));
cfg        = [];
cfg.trials = setdiff( data.trialinfo(:,end), tmp);
data       = ft_selectdata(cfg, data);
data.trialinfo    = data.trialinfo(:,1:(end-1));
data.cfg.previous = [];


%% save the identified artefacts
S.(var1) = toremove_overall0;
S.(var2) = toremove_muscle0;

if ~exist(art_filename) %no such a file at all
    save(art_filename, '-struct','S');
else
    save(art_filename, '-struct','S', '-append');
end

end