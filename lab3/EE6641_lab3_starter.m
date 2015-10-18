%   EE6641 ASAS Lab3: Spectrogram
%
%       starter.m
%
%   Created by Hao-Wei Lee, 2015.09

clear all; close all;

%windows
N=1000;
rectw=ones(N);
for i=1:N
hamm(i)=0.54-0.46*cos(2*pi*i/(N-1));
end
for i=1:N
han(i)=0.5*(1-cos(2*pi*i/(N-1)));
end

opt.useFunc = 1; % use built-in spectrogram function or not
opt.windowType = 'hann'; % rectwin, hann, or hamming etc.

D = 0.05;
overlap = 0.8;

if opt.useFunc, spectrogram = str2func('spectrogram');
else spectrogram = str2func('mySpectrogram'); end

window_gen = str2func( opt.windowType );

[sig, fs] = wavread('./test.wav');
win = window_gen( round(D*fs) ); %how many windows in a second

Nfft = power(2, ceil( log2(length(win)) )); %number of points(bins) in spectrogram
% Nfft = length(win);

S1 = spectrogram( sig(:, 1), win, round(length(win)*overlap), Nfft );
figure(1)
imshow( abs(S1(end:-1:1, :)) );

