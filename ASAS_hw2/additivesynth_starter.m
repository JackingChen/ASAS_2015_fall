function [y,finalState,temp_freqs] = additivesynth_starter(amps,freqs,N,initState,m,fs)
%%
% initState=state;

%%
J= size(amps,1); % # of tracks

finalState = zeros(J,3);
finalState(:,1) = amps(:,end);
finalState(:,2) = freqs(:,end);

y=zeros(N,1);

phi=zeros(J,N);
phi(:,1) = initState(:,3);
temp_amp = zeros(J,N);
temp_freqs = zeros(J,N);
beta= (1:N)/N;

%%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    temp_amp(:,i)=beta(i)*initState(:,1)+(1-beta(i))*finalState(:,1);
    temp_freqs(:,i)=beta(i)*initState(:,2)+(1-beta(i))*finalState(:,2);
end
temp_amp=10.^(temp_amp/20);
for i=2:N
    phi(:,i)=phi(:,i-1)+2*pi*temp_freqs(:,i)/fs;
end

for i=1:N
    y(i,:)=sum(temp_amp(:,i).*cos(2*pi*temp_freqs(:,i)+phi(:,i)));
end
%%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%


        


end