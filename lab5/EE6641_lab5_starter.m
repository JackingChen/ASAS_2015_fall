%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EE6641 Lab5                       %
% EE6641_lab5_starter.m             %
% created by Jeffrey Huang, 10/2015 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
close all;

%% Parameters
fs = 8000;
window_time = 0.05;
nBits = 16;
nChannels = 1;
poly_n=-1;
opt = 1;

%% Main
window_length = round(fs*window_time);

recorder = audiorecorder(fs,nBits,nChannels);
recObj = audiorecorder;

disp('Start recording.')
recordblocking(recObj, 2);
disp('End of Recording.');

fprintf('Playing your sound...\n');
play(recObj);

data = getaudiodata(recObj);
%%
if opt == 1
    result = spectrum_analysis(data,fs,window_length,poly_n);
else
    result = ACF(data,fs,window_length);
end

time_array = zeros(length(result.pitch),1);

for ii = 2:length(time_array)
    time_array(ii) = time_array(ii-1) + window_time/2; 
end

figure()
plot(time_array,result.pitch);
%%