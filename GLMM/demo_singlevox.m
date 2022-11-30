% Using the Pred_Alpha pipeline as an example.
% JY (Nov 2022)

% ========== Concatenate data =============
pIdx = 1; %parcel number 1, arbitrary for illustration purposes
tIdx = 1; %time point number 1, arbitrary
X    = [];
for iSub = 1:nsubjs
  f1 = sprintf( 'Sub%d_MEG', iSub); %behavioral data
  load( f1, 'results' );
  f2 = sprintf( 'Sub%d_sourcealpha_indvpeak_corticalsheet',iSub ); %alpha power
  load( f2, 'alphapar' );
  vec_yes      = (results.response == 1); %1 for target present, 0 for target absent
  vec_presence = results.presence;
  vec_alpha    = alphapar.trial(:,pIdx,tIdx); %alphapar.trial is rpt_chan_time in its dimension
  X = vertcat(X, horzcat( vec_yes, vec_presence, zscore(log(vec_alpha)), repmat(iSub, size(vec_alpha))) );
end

% ============= Apply GLMM =================
t = array2table(X, 'VariableNames',{'y', 'presence', 'alphapow', 'subject'});
    
% Test the effect of alphapow on y (i.e., probability report yes) while
% controling for stimulus presence and test for interaction effects.
% Model random subject intercepts and slopes.
% Intercept is modeled implicitly both at the group and subject level.
t.presence = categorical( t.presence ); %dummy-code
t.subject  = categorical( t.subject ); %dummy-code
m = fitglme( t, 'y ~ presence * alphapow + (presence*alphapow | subject)',...
        'Distribution', 'Binomial',...
        'Link', 'probit');
    
[PVAL,F,DF1,DF2] = coefTest( m, [0 0 0 1] );
fprintf( '\ndprime effect: pval = %s', num2str(PVAL) );
F_Dp = F;
    
[PVAL,F,DF1,DF2] = coefTest( m, [0 0 1 1/2] );
fprintf( '\ncriterion effect: pval = %s', num2str(PVAL) );
F_Crit = F;
