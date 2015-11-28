%% EE6641 HW: Linear prediction and Levinson-Durbin k-parameter estimation
% Created May 2013 as a homework.
% Last updated Nov 2015 for this year's Lab6 and HW3.
% Yi-Wen Liu
clear; close all;

DIR = './';
FILENAME = 'a.wav';

%FILENAME = 'EE6641-hw4-i.wav';
[y,fs1] = audioread([DIR FILENAME]);

soundsc(y,fs1);
fs = 16000;

y = resample(y,fs,fs1);
%% Parameters to play with
framelen = 0.05; % second. Please try changing this.
p = 20; % linear prediction order. Please try changing this.

%%
L = framelen*fs;

if L<=p
    disp('Linear prediction requires the num of equations to be greater than the number of variables.');
end

sw.emphasis = 1; % default = 1

numFrames = floor(length(y)/L);
excitat = zeros(size(y));
e_n = zeros(p+L,1);

LPcoeffs = zeros(p+1,numFrames);
Kcoeffs = zeros(p,numFrames); % reflection coeffs

Nfreqs = 1024; % Num points for plotting the inverse filter response
df = fs/2/Nfreqs;
ff = 0:df:fs/2-df;

if sw.emphasis == 1,
    y_emph = filter([1 -0.95],1,y);
else
    y_emph = y;
end

%% Linear prediction and estimation of the source e_n
win = ones(L,1); % Rectangular window.
for kk = 1:numFrames
    ind = (kk-1)*L+1:kk*L;
    ywin = y_emph(ind).*win';
    %A = lpc(ywin,p); %% This is actually the more direct way to obtain the
    % LP coefficients. But, in this script, we used levinson() instead
    % because it gives us the "reflection coefficients". 
    % 
    % (copied and modified from MATLAB's lpc() function)
    Y = fft(ywin,2^nextpow2(2*size(ywin,1)-1));
    R = ifft(abs(Y).^2);
    
    [A,errvar,K] = levinson(R,p);
    
    if kk == 1,
        e_n(p+1:end) = filter(A,[1],ywin);
    else
        ywin_extended = y((kk-1)*L+1-p:kk*L);
        e_n = filter(A,[1],ywin_extended);
    end
    excitat(ind) = e_n(p+1:end);
    figure(1);
    subplot(311);
    plot(ind/fs*1000, y(ind));
    xlabel('ms')
    set(gca,'xlim',[kk-1 kk]*framelen*1000);
    subplot(312);
    plot(ind/fs*1000, e_n(p+1:end));
    set(gca,'xlim',[kk-1 kk]*framelen*1000);
   
    subplot(313);
    [H,W] = freqz(1,A,Nfreqs);
    Hmag = 20*log10(abs(H));
    Ymag = 20*log10(abs(Y(1:Nfreqs)));
    Hmax = max(Hmag);
    offset = max(Hmag) - max(Ymag);
    plot(ff,Hmag); hold on;
    plot(ff,Ymag+offset,'r'); hold off;
    set(gca,'xlim',[0 fs/2],'ylim',[Hmax-50, Hmax+5]);
    xlabel('Hz')
    drawnow;
    %pause;
    LPcoeffs(:,kk) = A;
    Kcoeffs(:,kk) = K;
end
% play the estimated source signal
soundsc(excitat,fs); 
% show how cosistent the K-coefficients are.
[mean(Kcoeffs,2) sqrt(var(Kcoeffs,0,2))]

%% Below is a demo for vocal-tract shape estimation, part of HW3.
sw.ShapeEstim = 0;
if sw.ShapeEstim,
vtshape = zeros(p+1,numFrames); % defined as the sqrt of corss-sec area.
vtshape(end,:) = ones(1,numFrames); % Assume that the cross-section 
                                    % area near the glottis is fixed.
                                    
speed = 34300; % cm/s, speed of sound
dx = speed/fs; % cm per sample, spatial 
xx = (p-1)*dx:-dx:0;

for ii = 1:numFrames
    K = Kcoeffs(:,ii);
    for ll = p:-1:1
        refl = K(ll); % reflectance at this segment
        vtshape(ll,ii) = vtshape(ll+1,ii)*sqrt((1+refl)/(1-refl));
    end
    figure(4)
    plot(xx,vtshape(2:end,ii)); % First reflectance K(1) tends to be large because
            % impedance drastically changes between the end of vocal tract
            % and the open air outside.
    set(gca,'ydir','reverse');
    xlabel('glottis <-> lips (cm)')
    hold on;
    drawnow;
end
title(FILENAME)
end % if sw.ShapeEstim