function jy_batch_runICA( cfgcomp_filename, data )

cfg            = [];
cfg.channel    = 'MEG';
cfg.doscale    = 'no';
cfg.method     = 'runica';
cfg.runica.stop= 1e-8;
comp = ft_componentanalysis( cfg, data );

cfgcomp          = [];
cfgcomp.unmixing = comp.unmixing;
cfgcomp.topolabel= comp.topolabel;


if exist( cfgcomp_filename, 'file' )
    save( cfgcomp_filename, 'cfgcomp', '-append' );
else
    save( cfgcomp_filename, 'cfgcomp');
end


end