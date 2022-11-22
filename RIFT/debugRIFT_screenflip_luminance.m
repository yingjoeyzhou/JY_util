
clearvars; close all; clc;
addpath(genpath(pwd));

global devMode
devMode  = true; % false; %
freqRIFT = 1440;


durRIFT = 8; %in seconds


%% =========================================================== 
%% ===              Ask for user's inputs                  ===
%% =========================================================== 
[SubjectInfo, EnvInfo] = initExpt( 'PredOrLo', devMode, freqRIFT );




%% =========================================================== 
%% ===              Define the experiment                  ===
%% ===========================================================
global wPtr scr stim ctl


%% ============ non-MEG stuff ============
% === Open screen ===
[wPtr, wRect] = PsychImaging('OpenWindow', EnvInfo.screenID, [],...
                    EnvInfo.rectDisplay, [], [], [], 8);
Screen('BlendFunction', wPtr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% === The screen ===
scr              = [];
scr.rect         = wRect;
scr.idx          = EnvInfo.screenID;
scr.black        = BlackIndex(wPtr);
scr.white        = WhiteIndex(wPtr);
scr.gray         = (scr.white + scr.black)./2;
scr.ifi     = Screen('GetFlipInterval', wPtr);
scr.size    = EnvInfo.sizeDisplay;
scr.bgColor = repmat(scr.gray, [1,3]);
scr.viewDist= EnvInfo.viewDistance; %in cm
% double check IFI
if round(1 / scr.ifi) ~= 120
    warning('The IFI suggests a refresh rate of something else than 120 Hz!');
end
% the center pixel
[scr.xCtr, scr.yCtr] = RectCenter(scr.rect);
% the scale (pixels per degree)
cfg = []; cfg = scr; cfg.degrees = [1,1];
scr.ppdXY = JY_VisExptTools('deg2Pixel',cfg);
scr.ppdX  = mean(scr.ppdXY);
scr.ppdY  = mean(scr.ppdXY);
                
% === Define the stimulus (display duration and location) ===
stim = struct('Fix',[], 'dur',[], 'loc',[], 'patch',[]);

% Fixation
stim.Fix.color= repmat(scr.black, [1,3]);
stim.Fix.size = round(0.6 * scr.ppdX);
stim.Fix.type = 'bulleye'; %"cross", "dot" or "bulleye"
stim.Fix.mask = ones(stim.Fix.size+2, stim.Fix.size+2, 2); %+2 to ensure no obvious edge cause by cropping
stim.Fix.mask(:,:,1) = stim.Fix.mask(:,:,1) .* scr.gray;
stim.Fix.mask(:,:,2) = zeros( size(stim.Fix.mask(:,:,2)) );

stim.Fix.rectFix = CenterRectOnPoint([0 0 5 5], scr.xCtr, scr.yCtr);
     
% Location and size
stim.loc.xCtrDeg  = 3.8; %from Zhigalov & Jensen (2020) Hum Brain Mapp
stim.loc.xCtrPix  = stim.loc.xCtrDeg * scr.ppdX;
stim.loc.yCtrDeg  = 3.8; %from Zhigalov & Jensen (2020) Hum Brain Mapp
stim.loc.yCtrPix  = stim.loc.yCtrDeg * scr.ppdY;
stim.loc.rectDeg  = [0 0 5.7 5.7]; %from Zhigalov & Jensen (2020) Hum Brain Mapp
stim.loc.rectPix  = ceil( stim.loc.rectDeg .* scr.ppdX );
stim.loc.rect_R   = CenterRectOnPoint(stim.loc.rectPix, scr.xCtr + stim.loc.xCtrPix, stim.loc.yCtrPix + scr.yCtr);
stim.loc.rect_L   = CenterRectOnPoint(stim.loc.rectPix, scr.xCtr - stim.loc.xCtrPix, stim.loc.yCtrPix + scr.yCtr);

% define my bandpass-filtered patch
stim.patch.sizedeg   = stim.loc.rectDeg(3);
stim.patch.sizepix   = ceil( stim.patch.sizedeg .* scr.ppdX );
stim.patch.freq_mean = 2.5; %in cycles/deg
stim.patch.freq_sd   = 0.5; %0.5; %stim.patch.freq_mean / 3; %JY: arbitrary
stim.patch.ori_kappa = 400;
stim.patch.ori_mean  = []; %titrated
stim.patch.patchlum  = 0.5000;
stim.patch.patchcon  = 1; %full contrast
stim.patch.gauss_mask= 0;
stim.patch.gauss_sd  = [];
stim.noise           = stim.patch;

% define my aperture 
m.maskSiz = stim.patch.sizepix;
m.outerR  = stim.patch.sizedeg./2;
m.ppd     = scr.ppdX;
m.degSmoo = 0.5; %in dva
mm        = JY_VisExptTools('make_smooth_circular_mask', m);
stim.aperture  = 1-mm;


% === Define the experimental procedure ===
proc                 = [];
proc.nTrialsPerBlock = 64;
proc.nBlocks         = 8;
proc.nTrialsTotal    = proc.nTrialsPerBlock * proc.nBlocks;


%% ======= MEG-specific stuff ======
% === Define the rapid invisible frequency tagging (RIFT) ===
global RIFT 
RIFT           = [];
RIFT.broadband = 'yes';
RIFT.freq1     = 55;
RIFT.freq2     = 65;
RIFT.fs        = freqRIFT;
RIFT.totaldur  = durRIFT;
RIFT.timevec   = 0:(1/RIFT.fs):(RIFT.totaldur-1/RIFT.fs);
RIFT.fvec_S1L  = Shuffle( repmat(['1';'2'],[proc.nTrialsTotal./2, 1]) );
RIFT.fvec_S1R  = num2str( 3 - str2num( RIFT.fvec_S1L ) );

RIFT.Fix      = stim.Fix;
RIFT.Fix.size = ceil( RIFT.Fix.size ./ 2 );

RIFT.rectPix  = stim.loc.rectPix ./ 2;
RIFT.xCtrQuad = [ scr.xCtr.*0.5, scr.xCtr.*1.5, scr.xCtr.*0.5, scr.xCtr.*1.5 ];
RIFT.yCtrQuad = [ scr.yCtr.*0.5, scr.yCtr.*0.5, scr.yCtr.*1.5, scr.yCtr.*1.5 ];

for iQuad = 1:4 %becuase the the full screen has been divided into 4
    RIFT.rect_L{iQuad} = CenterRectOnPoint(RIFT.rectPix, RIFT.xCtrQuad(iQuad)-stim.loc.xCtrPix./2, RIFT.yCtrQuad(iQuad)+stim.loc.yCtrPix./2);
    RIFT.rect_R{iQuad} = CenterRectOnPoint(RIFT.rectPix, RIFT.xCtrQuad(iQuad)+stim.loc.xCtrPix./2, RIFT.yCtrQuad(iQuad)+stim.loc.yCtrPix./2);
end

% === After dividing the 1920x1080 frame into four 960x540 frames ===
% The number of pixels representing the stimulus is halfed
stim.patch.sizepix = ceil( stim.patch.sizepix ./ 2 );
stim.noise         = stim.patch;
% Same logic for the aperture, hence the same operation
m.maskSiz = stim.patch.sizepix;
m.outerR  = stim.patch.sizedeg./2;
m.ppd     = scr.ppdX ./ 2; %because of propixx cutting the screen into 4
m.degSmoo = 0.8; %in dva
mm        = JY_VisExptTools('make_smooth_circular_mask', m);
stim.outerpart = mm;
stim.aperture  = 1-mm;
% Place holders
stim.PH.img = ones( stim.patch.sizepix ) .* scr.gray;
stim.PH.img = stim.PH.img .* stim.aperture + scr.gray;


% ========= Define the patch for photodiode measurement =========
RIFT.patch           = [];
RIFT.patch.sizedeg   = 2;
RIFT.patch.sizepix   = ceil( RIFT.patch.sizedeg .* scr.ppdX ./ 2 );
RIFT.patch.xCtr      = [scr.xCtr, scr.rect(3), scr.xCtr, scr.rect(3)] - RIFT.patch.sizepix/2; %put it to the right bottom corner
RIFT.patch.yCtr      = [scr.yCtr, scr.yCtr, scr.rect(4), scr.rect(4)] - RIFT.patch.sizepix/2; %put it to the right bottom corner
RIFT.patch.rectPix   = [0 0 RIFT.patch.sizepix, RIFT.patch.sizepix];
RIFT.patch.locRects  = cell(4, 1);
for q = 1:4
    RIFT.patch.locRects{q} = CenterRectOnPoint(RIFT.patch.rectPix, ...
                                    RIFT.patch.xCtr(q), RIFT.patch.yCtr(q));
end
mR.maskSiz = RIFT.patch.sizepix;
mR.outerR  = RIFT.patch.sizedeg ./ 2;
mR.ppd     = scr.ppdX ./ 2; %because of propixx cutting the screen into 4
mR.degSmoo = 0.2; %in dva
mmR        = JY_VisExptTools('make_smooth_circular_mask', mR);
RIFT.patch.img  = ones( RIFT.patch.sizepix ) .* scr.gray;
RIFT.patch.aper = 1 - mmR;
RIFT.patch.img0 = RIFT.patch.img .* (1 - mmR) + scr.gray; %aperture = 1-mmR;






%% =========================================================== 
%% ===               photodiode placement                  ===
%% ===========================================================
Screen('FillRect', wPtr, [repmat(scr.gray,[1,3]),1], scr.rect );
texDio = Screen('MakeTexture', wPtr, RIFT.patch.img0);
texPH  = Screen('MakeTexture', wPtr, stim.PH.img);
for k = 1:4
    Screen('DrawTextures', wPtr, texDio, [], RIFT.patch.locRects{k});
    Screen('DrawTextures', wPtr, [texPH, texPH], [], [RIFT.rect_L{k}; RIFT.rect_R{k}]');
end
Screen('Flip', wPtr); pause(1); %KbWait();




%% =========================================================== 
%% ===              show stimuli w/ RIFT                   ===
%% ===========================================================
trialcfg = [];
trialcfg.freqL = RIFT.freq1; %JY: arbitrary
trialcfg.freqR = RIFT.freq2; %JY: arbitrary

trialcfg.S1Ori = 'V'; %proc.S1Ori(nCount);
trialcfg.S2Ori = 'H'; %proc.S2Ori(nCount);

trialcfg.S1 = genBandpassGratingRaw( stim.patch, trialcfg.S1Ori ); %JY: we will put this on the left
trialcfg.S2 = genBandpassGratingRaw( stim.patch, trialcfg.S2Ori ); %JY: we will put this on the right

switch RIFT.broadband
    case 'yes'
        trialcfg.S1TagL = defineBroadbandLuminance( RIFT.timevec, trialcfg.freqL );
        trialcfg.S1TagR = defineBroadbandLuminance( RIFT.timevec, trialcfg.freqR );
    case 'no'
        trialcfg.phaseL = rand().*pi;
        trialcfg.phaseR = trialcfg.phaseL + pi; %anti-phase
        trialcfg.S1TagL = sin(2*pi .* trialcfg.freqL .* RIFT.timevec + trialcfg.phaseL);
        trialcfg.S1TagR = sin(2*pi .* trialcfg.freqR .* RIFT.timevec + trialcfg.phaseR);
end


%% inside the trialfun: v1
%{ %
tmp = stim.aperture;
nFlips = durRIFT .* 120; %JYL hard-coded 120 Hz, which is the graphic card's refresh rate

% reshape the freq-tagged signal (i.e., contrast)
lumL = reshape( trialcfg.S1TagL, [12, numel(trialcfg.S1TagL)./12] ) .*0.5 + 0.25;
lumR = reshape( trialcfg.S1TagR, [12, numel(trialcfg.S1TagR)./12] ) .*0.5 + 0.25;

trialcfg.S1.img = trialcfg.S1.img .* 0.5;% .* stim.aperture;% + stim.patch.patchlum;
trialcfg.S2.img = trialcfg.S2.img .* 0.5;% .* stim.aperture;% + stim.patch.patchlum;


tmp = [];
for iFlip = 1:nFlips
    
    tmpL = lumL(:,iFlip);
    tmpR = lumR(:,iFlip); 
    tmpL = reshape( tmpL, [1, 4, 3] ); %from 1st to 3rd column, Reds -> Greens -> Blues;
    tmpR = reshape( tmpR, [1, 4, 3] ); %from 1st to 3rd column, Reds -> Greens -> Blues;
    
    for k = 1:4 %The four quardrants
        % imgL = bsxfun( @times, trialcfg.S1.img, tmpL(1,k,:) ) + stim.patch.patchlum;
        % imgR = bsxfun( @times, trialcfg.S2.img, tmpR(1,k,:) ) + stim.patch.patchlum;
        % dioR = bsxfun( @times, RIFT.patch.img .* RIFT.patch.aper, tmpR(1,k,:) ) + stim.patch.patchlum;
        imgL = (trialcfg.S1.img + tmpL(1,k,:)).* stim.aperture + stim.outerpart.*stim.patch.patchlum;
        imgR = (trialcfg.S2.img + tmpR(1,k,:)).* stim.aperture + stim.outerpart.*stim.patch.patchlum;
        dioR = bsxfun( @times, RIFT.patch.img .* RIFT.patch.aper, tmpR(1,k,:) ) + stim.patch.patchlum;
        texL = Screen('MakeTexture',wPtr,imgL);
        texR = Screen('MakeTexture',wPtr,imgR);
        txDioR = Screen('MakeTexture',wPtr,dioR);
        Screen('DrawTextures',wPtr,[texL,texR,txDioR],[],...
            [RIFT.rect_L{k}; RIFT.rect_R{k}; RIFT.patch.locRects{k}]',[],[]);
    end
    
    vbl = Screen('Flip',wPtr);
    tmp = [tmp, vbl];
    
    if devMode, pause(1); end
    
    if iFlip==nFlips; break; end
    
end
%}



sca;
Datapixx('SetPropixxDlpSequenceProgram', 0); %Revert to standard 120Hz refresh rate
Datapixx('RegWrRd');
Datapixx('Close');

save( SubjectInfo.datafile );