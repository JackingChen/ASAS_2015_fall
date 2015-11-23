clear; close all;
sw.plot = 0;
% [x,fs] = wavread('wrenpn.wav');
% [x,fs] = wavread('peaches_16.wav');
[x,fs] = wavread('mymom_16.wav');
% [x,fs] = wavread('draw_16.wav');
x = x/max(abs(x)); sound(x,fs);
nx = length(x);

%% ANALYSIS PARAMETERS
frameRate =40;  % frames per second
M = floor(fs/frameRate);  
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%
%w = hanning(M);
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

for m=1:nFrames-1
%     %%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
    [amp,freq] = findpeaks_starter_sy(spectra{m}(:,2),spectra{m}(:,1),maxPeaks);
    amps(:,m) = amp;
    freqs(:,m) = freq;
    %%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
end
aaa=freqs;
bbb=amps;
% figure(1)
% plot(8000/pi*freqs');
% xlabel('frame');
% ylabel('hz');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTHESIS
R = round(R* expandRatio);  % time expansion
% freqs = min(pi,freqs*fRatio);
y = zeros((nFrames+1)*R,1);
state = zeros(maxPeaks,3);  % [ampInitials, freqInitials, phaseInitials] 
% state(:,1) = ... ;
state(:,2) = freqs(:,m);

for m=1:nFrames-1
    %%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
    rr=R*(m-1)+1:1:R*(m-1)+R;
    y_synth = additivesynth_starter_sy(amps(:,m),freqs(:,m),R,state,rr);
    y(rr') = y_synth;
    %%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
end

y = y/(max(abs(y))+1); 
sound(y,fs);
% wavwrite(y,fs,'LinearMethod.wav');



