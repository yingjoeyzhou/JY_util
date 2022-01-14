function data = jy_checkcomponent( preprocLogFile, data)        

% load the cfgcomp (output of jy_batch_runICA)
load(preprocLogFile,'cfgcomp');


% project the components back to the original data
cfg           = [];
cfg.doscale   = 'no';
cfg.unmixing  = cfgcomp.unmixing;
cfg.topolabel = cfgcomp.topolabel;
comp = ft_componentanalysis(cfg, data);


% visualize the components
cfg = [];
cfg.layout   = 'CTF275.lay';
cfg.viewmode = 'component';
cfg.zlim     = 'maxabs';
ft_databrowser(cfg, comp);


% indicate the components to remove
rmcomp = str2num(input('\n\n\n ***** Indicate the components to remove: ****** \n','s'));


% save the to-be-removed components
save(preprocLogFile, 'rmcomp', '-append');


% reject the components
cfg = []; cfg.component = rmcomp;
data= ft_rejectcomponent(cfg, comp, data);

end