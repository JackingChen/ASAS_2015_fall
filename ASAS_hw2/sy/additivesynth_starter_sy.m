function y = additivesynth_starter_sy(amps,freqs,N,initState,R)

J= length(amps); % # of tracks
fs=16000;
finalState = zeros(J,3);
% finalState(:,1) = amps;
% finalState(:,2) = freqs;
y=zeros(1,N);
phi = initState(:,3);
temp_amp = zeros(1,N);
temp_freqs = zeros(1,N);
w = blackman(N);
% sinn=zeros(4000,1);

%%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
for jj =1:J
    y=y+amps(jj).*w'.*sin(2*pi*freqs(jj).*R/fs); %還要再多一個window function
%%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
end
y=y';
end