%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EE6641 HW1                        %
% sound_generator.m                 %
% created by Jeffrey Huang, 10/2015 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
close all;

%% Parameters
opt = 0;
Dt=3;
Df=1/(4*pi)/Dt;

dt = 0.1;
df = 1/(4*pi)/dt;
f1 = 110;
f2 = 450;

%% tasks
[correctness ucert] = sound_generator(Dt,Df,dt,df,f1,f2,opt);

%% Plotting
if isstruct(correctness)
    test_times = correctness.stage_idx;
    figure()
    scatter(1:test_times,correctness.time(1:test_times));
    title('Your hearing resolution in time')
    xlabel('(times)');
    ylabel('(sec) (negative means incorrect)')
    
    figure()
    scatter(1:test_times,correctness.frequency(1:test_times));
    title('Your hearing resolution in frequency')
    xlabel('(times)');
    ylabel('(hz) (negative means incorrect)')
    fprintf('Your uncertainty limit is %.7f \n',4*pi*correctness.resolution(1)*correctness.resolution(2));
end