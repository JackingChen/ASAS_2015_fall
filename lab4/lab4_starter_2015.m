% EE6641 Analysis and Synthesis of Audio Signals 
%
% Homework #2 Problem 2 starter code
%
% FIR filter design by windowing, and sequential processing using the 
% Overlap-Add (OLA) analysis and synthesis framework
% 
% Created Oct. 19, 2010. Updated March 11, 2013.
% Sencod edition in March 19 2014
%  

clear; close all;
sw.plotFFT = 0; %Plot the spectrum if sw.plotFFT = 1

DIR = './';
[x, fs] = wavread([DIR 'boxing.wav']);
    % ** STUDENTS: find your favorite song and convert it to a .wav file
x = x(1:10*fs); x = x(:); % take the first 10 seconds
L = length(x);
N = 512;
hopsize = N/2;
win = hann(N+1); % Hann window
win = win(1:end-1);

num_frame = 2*floor(L/N)-2;
x = x(1:(num_frame+1)*(N/2));
x = x(:);


% THE code below designs a low-pass filter using fir1().
%
% *** STUDENTS: DESIGN A BAND-STOP FILTER INSTEAD.  The goal is to remove
% human voice while keeping as much of background music as possible.
% Recommended cut-off frequencies are 200 Hz and 2500 Hz, but the actual
% effectiveness of voice-removal may depend on the pitch range of the
% singer (male/female/children), type of background music, etc.
%
% type 'help fir1' for more details.
% 
f_cut_head = 2200; % Hz
f_cut_head_norm=double(f_cut_head/(fs/2));
f_cut_tail=  2500; 
f_cut_tail_norm= double(f_cut_tail/(fs/2));
h = fir1(128,[f_cut_head_norm f_cut_tail_norm]);
h = h(:);

% set length of output
N_taps = length(h);
y = zeros(length(x)+N_taps-1,1); %Length of the ouput signal = N+P-1

% Set FFT number of bins
N_zp = 1024;

if sw.plotFFT,
    df = fs/N_zp;
    ff = 0:df:(N_zp-1)*df;   
end

% Reserve memory space
x_zp = zeros(N_zp,1);
h_zp = zeros(N_zp,1);
h_zp(1:N_taps) = h; % zero-padded impulse response
H = fft(h_zp);
%%
%********** Please implement your OLA algorithm *************%
tic
y2=[];
for mm = 1:num_frame
    t_start = (mm-1)*hopsize;    
    tt = (t_start+1):(t_start+N);
    x_win = x(tt).*win; % windowing
    x_zp= zeros(N_zp,1); % zero-padding
    x_zp(1:length(x_win))=x_win;
    
    % Filtering using direct multiplication in freq. domain
    X = fft(x_zp); 
    Y = H.*X; 
    y_win = ifft(Y); %ifft
    
    % plot FFT
    if sw.plotFFT,
        figure(1);
        Xmag = 20*log10(abs(X));
        Ymag = 20*log10(abs(Y));
        plot(ff/1000,Xmag,'b',ff/1000,Ymag,'r--');
        set(gca,'xlim',[0 fs/2/1000],'ylim',[-100 20]);
        
        xlabel('kHz');
        ylabel('dB');
        title(sprintf('t = %.3f s\n',t_start/fs));
        drawnow;
    end
    
    % overlap-add 
    tt2 = 1:(N+N_taps-1);
    y(t_start+tt2) = y(t_start+tt2)+y_win(tt2); 
    %y2=[y2; y_win(tt2)];
    
    % plot FFT
    if sw.plotFFT,
        figure(1);
        Xmag = 20*log10(abs(X));
        Ymag = 20*log10(abs(Y));
        plot(ff/1000,Xmag,'b',ff/1000,Ymag,'r--');
        set(gca,'xlim',[0 fs/2/1000],'ylim',[-100 20]);
        
        xlabel('kHz');
        ylabel('dB');
        title(sprintf('t = %.3f s\n',t_start/fs));
        drawnow;
    end
end
toc;
%************************************************************%
%sound(y,fs);
wavwrite(y,fs,'lab4.wav');
