%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EE6641 HW2                           %
% Time-warping and frequency shifting  %
% based on sinusoidal modeling.m       %
% created by Jeffrey Huang, 11/2015    % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
sw.plot = 0;
% [x,fs] = wavread('wrenpn.wav');
% [x,fs] = wavread('peaches_16.wav');
[x,fs] = wavread('mymom_16.wav');
% [x,fs] = wavread('draw_16.wav');
x = x/max(abs(x)); sound(x,fs);
nx = length(x);



%% ANALYSIS PARAMETERS
frameRate =40;  % frames per second
M = floor(fs/frameRate); % tmp parameter
nFrames = floor(nx/M)*2-1;
R = floor(M/2);  % Exact COLA not required
N = 2^(1+floor(log2(5*M+1))); % FFT length, at least a factor of 5 zero-padding

maxPeaks = input('how many peaks to track? '); 
	% Up to this many sinusoidal peaks (birdsong-specific!)
expandRatio = input('time expansion factor?');
freqShift = input('how many semitones of pitch shift?');
fRatio = 2^(freqShift/12); 

%% VECTOR VARIABLES DECLARATION
amps = zeros(maxPeaks,nFrames);
freqs = zeros(maxPeaks,nFrames);
% x_f=zeros(1,frame_size);
x_zp=zeros(1,N);
w_zp=zeros(1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%
%w = hanning(M);

df = fs/N;
ff = 0:df:(N-1)*df;
frame_size=floor(nx/nFrames);
% hopsize=180;
w = blackman(frame_size);
spec=zeros(frame_size,nFrames);

% demonstrate%%%%%%%%%%%%%%%%%%%%%%%%
% s=specgram(x,512,fs);
% [R,M ]= extractmax(abs(s));
% disp(['size of R is ',num2str(size(R,1)),'rows x',numm2str(R,2),' cols']);
% tt= [1:size(R,2)]*128/fs;
% F=R*fs/256;
% specgram(x,256,fs);
% colormap(1-gray);
% hold on
% plot(tt,F','r');
threshold=0.1;
spec=spectrogram(x,w,frame_size/2,512,fs,'yaxis');
spec=abs(spec);
for m=1:nFrames
    %%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
%     x_zp(1:frame_size)=x((m-1)*frame_size+1:m*frame_size);
%     x_f=fft(x_zp);
%     w_f=fft(w);
%     w_zp(1:frame_size)=w_f;
%     x_f=x_f.*w_zp;
%     x_f=x_f(1:frame_size);
%     x_f=abs(x_f);
    [amps(:,m),freqs(:,m)]=findpeaks_starter(spec(:,m),maxPeaks,threshold);

%     x_f=abs(x_f);
%     spec(:,m)=x_f;
    %%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
end


aaa=freqs;
bbb=amps;
% spec=spectrogram(x,w,0,512,fs,'yaxis');
% spec=abs(spec);
freqs=freqs*(pi/8000);

% [S,F,T]=spectrogram(x,w,0,512,fs);
% S=abs(S);
% imagesc(S/max(abs(S(:))) );
figure
spectrogram(x,w,0,512,fs,'yaxis');
set(gca,'YDir', 'normal');
colormap(1-gray);
hold on;
xx=linspace(0,1,length(freqs));
% freqs=freqs;
plot(xx,8000/pi*freqs');

%
xlabel('frame');
ylabel('hz');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTHESIS
R = round(R* expandRatio);  % time expansion
freqs = min(pi,freqs*fRatio);
% freqs=freqs*8000/pi;
y = zeros((nFrames+1)*R,1);
state = zeros(maxPeaks,3);  % [ampInitials, freqInitials, phaseInitials] 
state(:,1) = amps(:,1) ;
state(:,2) = freqs(:,1) ;

line=[];
freqs_c=[];
for m=1:nFrames-1
%%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
    [y,finalState,temp_freqs] = additivesynth_starter(amps(:,m),freqs(:,m),frame_size,state,m,fs)
    line=[line; y];
    freqs_c=[freqs_c temp_freqs];
%%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
end
y=line;
figure
specgram(x,512,fs);
colorbar off ;colormap bone;
hold on;
xx=1:1:length(freqs_c(1,:));
plot(xx,freqs_c(1,:));
figure
plot(freqs);


y = y/(max(abs(y))+1); 
sound(y,fs);
wavwrite(y,fs,'LinearMethod.wav');



