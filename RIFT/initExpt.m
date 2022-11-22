function [SubjectInfo, EnvInfo] = initExpt( exptName, devMode, propixxFreq )
% Get subject & environment information by asking for user's inputs.
%   INPUT:
%       exptName   : string, name of the experiment.
%       devMode    : logical, true or false.
%       propixxFreq: scalar, refresh rate of the propixx projector.
% 
%   OUTPUT:
%       SubjectInfo: structure, fields contain info about the subject.
%       EnvInfo    : structure, fields contain info about the environment.
% 
% JY (Sep, 2022)

t = clock;

SubjectInfo            = [];
SubjectInfo.SubjectID  = input('Subject ID (PXX, P = version):\n','s');
SubjectInfo.Age        = str2double( input('Age:\n','s') );
SubjectInfo.Gender     = input('Gender (F/M/N, N for non-binary):\n','s');
SubjectInfo.Handedness = input('Handedness (L/R):\n','s');

SubjectInfo.timestamp = sprintf('%d_%d_%d_%d_%d_%d',t(1),t(2),t(3),t(4),t(5),round(t(6)));
SubjectInfo.datafile  = [exptName,'_Sub',SubjectInfo.SubjectID,'_',SubjectInfo.timestamp,'.mat'];
SubjectInfo.randseed  = t(4)*100 + t(5);

rng(SubjectInfo.randseed);


% Here we call some default settings for setting up Psychtoolbox
% A 'featureLevel' of 2 will execute the AssertOpenGL command, the
% KbName('UnifyKeyNames') command, and the
% Screen('CollorRange',window,1,[],1) command to switch the default color
% range from the 0-255 integer number to normalized floating point number
% range 0.0 - 1.0 to unify color specifications across displays/devices.
PsychDefaultSetup(2);
% Get the screen numbers
screens = Screen('Screens');
% Draw to the external screen if avaliable
screenNumber = max(screens);


% Output EnvInfo
EnvInfo              = [];
EnvInfo.screenID     = screenNumber;
EnvInfo.skipSyncTest = devMode;

if devMode %in developer mode, e.g., while debugging
    
   EnvInfo.sizeDisplay  = [52.7, 29.6]; %in centimeters
   EnvInfo.viewDistance = 57; %in centimeters
   EnvInfo.rectDisplay  = [0 0 960 540];
   EnvInfo.calibFile    = [];
   EnvInfo.propixx      = false;
   EnvInfo.eyelink      = false;
   
else %NOT in developer mode, i.e., running for real in MEG w/ propixx
    
   EnvInfo.sizeDisplay  = [56, 31.5]; %in centimeters
   EnvInfo.viewDistance = 115; %in centimeters
   EnvInfo.rectDisplay  = [];
   EnvInfo.calibFile    = [];
   EnvInfo.propixx      = true;
   EnvInfo.eyelink      = true;
   
   % ------ initialize DATAPixx here ---------
   isConnected = Datapixx('isReady');
   if ~isConnected
       Datapixx('Open');
       Datapixx('StopAllSchedules');
   end
   
   if isempty( propixxFreq ) | propixxFreq==1440
       %Set Propixx to 1440Hz refresh (also known as Quad12x)
       Datapixx('SetPropixxDlpSequenceProgram', 5);
       disp('RIFTing at 1440 Hz!');
       
   else
       Datapixx('SetPropixxDlpSequenceProgram', 2);
       disp('RIFTing at 480 Hz!');
   end
   
   % Synchronize DATAPixx registers to local register cache
   Datapixx('RegWrRd');
   
    
% %    % Show how many TTL output bits are in the Datapixx
% %    nBits = Datapixx('GetDoutNumBits');
% %    fprintf('\nDATAPixx has %d TTL output bits\n\n', nBits);
% %     
% %    % Bring all the outputs low
% %    Datapixx('SetDoutValues', 0);
% %    Datapixx('RegWrRd');
% %    disp('datapixx set up'); 
   
end

end