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
[x,fs] = wavread('draw_16.wav');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%
%w = hanning(M);

% df = fs/N;
% ff = 0:df:(N-1)*df;
% frame_size=floor(nx/nFrames);
% % hopsize=180;
% w = blackman(frame_size);
% spec=zeros(frame_size,nFrames);

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

w = blackman(M);
df = fs/N;
ff = 0:df:(N-1)*df;

window_length = M;
wav_file = x;
hop_size = window_length/2;
wav_length = length(wav_file);
st_num = nFrames;
ed_num = st_num;
frame_num = nFrames;
spectra = cell(frame_num,1);
st = zeros(st_num,1);
ed = zeros(ed_num,1);
window_vec = w;

for ii = 1:frame_num
    if ii == 1
        st(ii) = 1 + (window_length)*(ii-1);
        ed(ii) = window_length*ii;
    elseif ii == frame_num
        st(ii) = st(ii-1) + hop_size;
        ed(ii) = wav_length;
    else
        st(ii) = st(ii-1) + hop_size;
        ed(ii) = ed(ii-1) + hop_size;
    end
end

for ii = 1:frame_num
    if ii ~= frame_num                
        temp = abs(fft(wav_file(st(ii):ed(ii)).*window_vec));
        spectra{ii} = zeros(floor(floor(length(temp)))/2+1,2);
        spectra{ii}(:,1) = 0:fs/length(temp):fs/2;
        spectra{ii}(:,2) = temp(1:floor(length(temp))/2+1);
    else
        spectra{ii} = zeros(floor(length(temp))/2+1,2);
        spectra{ii}(:,1) = 0:(fs/length(temp)):(fs/2);
    end
end

for m=1:nFrames
    %%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
%     x_zp(1:frame_size)=x((m-1)*frame_size+1:m*frame_size);
%     x_f=fft(x_zp);
%     w_f=fft(w);
%     w_zp(1:frame_size)=w_f;
%     x_f=x_f.*w_zp;
%     x_f=x_f(1:frame_size);
%     x_f=abs(x_f);
    [amps(:,m),freqs(:,m)]=findpeaks_starter(spectra{m}(:,1),spectra{m}(:,2),maxPeaks);

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
figure(1)
spectrogram(x,w,0,512,fs,'yaxis');
set(gca,'YDir', 'normal');
colormap(1-gray);
hold on;
xx=linspace(0,1,length(freqs));
% freqs=freqs;
plot(xx,freqs');

%
xlabel('frame');
ylabel('hz');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTHESIS
R = round(R* expandRatio);  % time expansion
freqs = min(pi,freqs*fRatio);
freqs=freqs*8000/pi;
y = zeros((nFrames+1)*R,1);
state = zeros(maxPeaks,3);  % [ampInitials, freqInitials, phaseInitials] 
state(:,1) = amps(:,1) ;
state(:,2) = freqs(:,1) ;

line=[];
freqs_c=[];
for m=1:nFrames-1
%%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
    rr=R*(m-1)+1:1:R*(m-1)+R;
    [y_frame,finalState,temp_freqs] = additivesynth_starter(amps(:,m),freqs(:,m),R,state,fs,rr);
    y(rr')=y_frame;
    freqs_c=[freqs_c temp_freqs];
%%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
end

figure(2)
specgram(x,512,fs);
colorbar off ;colormap bone;
hold on;
xx=1:1:length(freqs_c(1,:));
plot(xx,freqs_c(1,:));
figure(3)
plot(freqs);


y = y/(max(abs(y))+1); 
sound(y,fs);
wavwrite(y,fs,'LinearMethod.wav');



