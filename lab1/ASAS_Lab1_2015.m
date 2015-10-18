% EE6641 Lab 1: Tone synthesis by directly creating the discrete-time 
% waveform
% 
% Sep 15, 2015
% Yi-Wen Liu

clear; close all
fs = 8000;
dt = 1/fs;
dur = 2;
f1 = [941,852,697,697,770,770,852,852,697,852]; % Hz
f2 = [1336,1477,1477,1209,1477,1477,1477,1336,1209,1336]; % Hz
A1 = 0.5;
A2 = 0.5;

numsamples = dur/dt;

tt = 0:dt:(numsamples-1)*dt;
tt = tt(:);%transpose

% for i =1:length(f1)
%     x = A1*cos(2*pi*f1(i)*tt) + A2*cos(2*pi*f2(i)*tt);
%     sound(x(1:3000),fs);
%     pause(0.2)
% 
%     
% end
finalS=[]
for i =1:10
    phi=i*pi/2/6
    x = A1*cos(2*pi*f1(i)*tt+phi) + A2*cos(2*pi*f2(i)*tt+phi);
    
    finalS=cat(1,finalS,x(1:4000));
    sound(x(1:9000),fs);
    pause(0.5)

    
end
audiowrite('~/Desktop/Hw0.wav',finalS,fs)
% figure(i)
% plot(tt*1000,x,'--');
% set(gca,'xlim',[0 5]);
% xlabel('msec');
% hold on;
% stem(tt(1:1:160)*1000,x(1:1:160))
% setFontSizeForAll(14);
x=0;

