%
% Compare the photodiode's recorded signal with the RIFTing signal used by
% PsychToolbox.
% 
% During the pilot, the photodiode was placed on top of a circle at the 
% bottom-right of the screen. We sent triggers to mark the onset of the 
% stimulus (code = 4). These codes should be used to epoch the data in the 
% .fif file. 
% 
% Each .mat file in this folder corresponds to a few 5s trials. We need to
% use the timestamps of the .mat files to figure out which trials were
% shown earlier and which were shown later.
% 
% JY (Nov 2022)

clearvars; close all; clc;
addpath( '/Users/joeyzhou/Documents/GitHub/fieldtrip' );
ft_defaults;


%% Epoch the .fif file
tit = '5k';
fif = dir(['*',tit,'_raw.fif']);
hdr = ft_read_header( fif.name );
dat = ft_read_data( fif.name );

% ====== Define the channels of interest ======
iDiodeA = cellfun(@(x) contains(x,'MISC011'), hdr.label, 'uniformoutput',true);
iDiodeB = cellfun(@(x) contains(x,'MISC012'), hdr.label, 'uniformoutput',true);
iTrigg  = cellfun(@(x) contains(x,'STI101'), hdr.label, 'uniformoutput',true);

% ====== The Nth sample where we have the trigger ===
idxS1  = strfind( dat(iTrigg,:), [0 4 4] ) + 1;

% ========= Trial length (in seconds) ==========
dur = 5;


%% Compare the signals via visualization
tFIF = (1/hdr.Fs):(1/hdr.Fs):dur;
segFIF = [1:(hdr.Fs .* 0.25)];% + hdr.Fs.*9.75;

yA = dat(iDiodeA, :);
yA = ( yA-min(yA) )./ range(yA);

yB = dat(iDiodeB, :);
yB = ( yB-min(yB) )./ range(yB);

for ii = 1:numel( idxS1 )
    iBeg = idxS1(ii);
    iEnd = idxS1(ii) + hdr.Fs*dur - 1;
    
    yDiodeA = yA(iBeg:iEnd);
    yDiodeB = yB(iBeg:iEnd);
    
    figure(1), clf, hold on,
    set(gcf,'units','centimeters', 'position', [0 0 50 15]);
    l1=plot( tFIF(segFIF), yDiodeA(segFIF), 'k.-', 'linewidth',1 );
    l2=plot( tFIF(segFIF), yDiodeB(segFIF), 'r.-', 'linewidth',1 );
    legend([l1,l2],{'photodiode S', 'photodiode L'}); legend boxoff;
    set(gca,'tickdir','out','fontsize',18);
    ylabel('signal (a.u.)', 'fontweight','bold');
    xlabel('time (s) from stimulus onset', 'fontweight','bold');
    ylim( [0, 1.2]);
    pause(3);
end


%% plot fourier spectrum
for ii = 1:numel(idxS1)
    iBeg = idxS1(ii);
    iEnd = idxS1(ii) + hdr.Fs*dur - 1;
    
    yDiodeA = yA(iBeg:iEnd);
    yDiodeB = yB(iBeg:iEnd);
    
    figure(1), clf, hold on,
    plotFFT( tFIF, yDiodeA, 'k' );
    plotFFT( tFIF, yDiodeB, 'r' );
    
    pause(3);
end



%% compute instantaneous phase
tPTB = (1/1440):(1/1440):dur;
tFIF = (1/hdr.Fs):(1/hdr.Fs):dur;

segPTB = [1:(1440 .* 0.25)];% + 1440.*9.75;
segFIF = [1:(hdr.Fs .* 0.25)];% + hdr.Fs.*9.75;

for ii = 1:numel( idxS1 )
    iBeg = idxS1(ii);
    iEnd = idxS1(ii) + hdr.Fs*dur - 1;
    
    yDiodeA = yA(iBeg:iEnd);
    yDiodeB = yB(iBeg:iEnd);
    
    zDiode = hilbert( yDiodeA );
    inst_phase_diodeA = angle(zDiode);%unwrap( angle(zDiode) );
    
    zDiode = hilbert( yDiodeB );
    inst_phase_diodeB = angle(zDiode);%unwrap( angle(zDiode) );
    
    figure(1), clf, hold on,
    set(gcf,'units','centimeters', 'position', [0 0 50 15]);
    l1 = plot( tFIF(segFIF), inst_phase_diodeA(segFIF), 'k.-', 'linewidth',1 );
    l2 = plot( tFIF(segFIF), inst_phase_diodeB(segFIF), 'r.-', 'linewidth',1 );
    legend([l1,l2], {'photodiode S', 'photodiodeL'}, 'orientation','horizontal'); legend boxoff;
    set(gca,'tickdir','out','fontsize',18);
    title( ['MEG sampling rate = ', tit] );
    ylabel('instantaneous phase (rad)', 'fontweight','bold');
    xlabel('time (s) from stimulus onset', 'fontweight','bold');
    ylim( [-1,1].* max(abs(inst_phase_diodeA(segFIF))) );
    
    
    figure(2), clf, hold on,
    plot( inst_phase_diodeA(segFIF), inst_phase_diodeB(segFIF), 'o');
    plot( [-1, 1], [-1, 1], 'k:','linewidth',2);
    tx = text( 0.5, -0.7, ['r = ',num2str( round(corr( inst_phase_diodeA', inst_phase_diodeB'),3) )] );
    tx.FontSize = 16;
    xlabel('Diode A', 'fontweight','bold');
    ylabel('Diode B', 'fontweight','bold');
    xlim( [-1,1] );
    ylim( [-1,1] );
    set(gca,'TickDir','out');
    pause;
    
end


%% ========= sub-function ============
function plotFFT( t, x, lnspec )

Fs= 1./diff(t(1:2));
t = t - Fs;
L = numel(t);


Y = fft(x);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

f2plot = [0, 100];
idxBeg = 1;
idxEnd = find( f<= f2plot(2), 1, 'last' );

plot(f(idxBeg:idxEnd),P1(idxBeg:idxEnd), lnspec) 
% title("Single-Sided Amplitude Spectrum of X(t)")
% xlabel("f (Hz)")
% ylabel("|P1(f)|")

end
