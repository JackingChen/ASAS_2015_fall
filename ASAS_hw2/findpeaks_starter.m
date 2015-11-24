function [amps,freqs]=findpeaks_starter(F,X,maxNumPeaks);
% % %%
% F=spectra{1}(:,1);
% X=spectra{1}(:,2);
% maxNumPeaks=maxPeaks;
%%
if max(X)==0
    amps=zeros(maxNumPeaks,1);
    freqs=zeros(maxNumPeaks,1);
    return
end
N=size(X(:),1);
data = 20*log10(abs(X)); %intensity in db scale

data=data(:);
ind=find( (data>[data(1)-100;data((1:N-1)')]) ...
    & (data>=[data((2:N)'); data(N)-100]) );   % find peaks. "ind" means the locations of the peaks

Peaksmax=length(ind);
% figure
% datapeak=zeros(1,length(data));
%     datapeak(ind)=data(ind);
%     xx=1:1:length(data);
%     plot(xx,data',xx,datapeak,'o');
%%
%
% You might need to use SORTROWS() to identify highest peaks for each
% spectrum
%%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
    [a ind_S]=sort(data(ind),'descend');
    
    
    Data=data(ind);
    Data_sort=sortrows(Data);
    Data_sort = Data_sort(end:-1:1);
    peaks_collect=[];
    if isempty(Data_sort)
        Data_sort=zeros(length(data),1);
    end
    % show data peak funciton
    
    
%%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
%%

if Peaksmax >= maxNumPeaks
    peaktmp=Data_sort(1:maxNumPeaks);
    amps=zeros(maxNumPeaks,1);
    freqs=zeros(maxNumPeaks,1);
    for i=1:maxNumPeaks
        indx=find(data==peaktmp(i));
        amps(i)=data(indx);
        freqs(i)=F(indx);
    end
    
else
    peaktmp=Data_sort(1:Peaksmax);
    Peak_lack = maxNumPeaks-Peaksmax; % Peak 缺少的個數
    delta=zeros(Peaksmax,1);
    
    for k=1:Peaksmax-1
        delta(k)=Data_sort(k+1)-Data_sort(k);
    end
    z=sort(delta,'descend');
    inter=zeros(Peak_lack,1);
    for bb=1:Peak_lack
        vv=find(delta==z(bb));
        inter(bb)=(Data_sort(vv)-Data_sort(vv+1))/2;
    end
    amps=zeros(maxNumPeaks,1);
    freqs=zeros(maxNumPeaks,1);
    for j=1:Peaksmax
            %% Implement quadratic interpolation here, and save the
        %%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
        indx=find(data==peaktmp(j));
        amps(j)=data(indx);
        freqs(j)=F(indx);
        %%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%spectra{30}(:,2)
    end
    YI=interp1(data,F,inter);
    freqs(Peaksmax+1:maxNumPeaks)=YI;
    amps(Peaksmax+1:maxNumPeaks)=inter;
% for j=1:maxNumPeaks
%         %% Implement quadratic interpolation here, and save the
%         % amplitude and frequencies of the peaks in a matrix peaks.
%         % The first column of PEAKS should be all the estimated amplitudes
%         % of sinusoidal components, and the second column of PEAKS should
%         % be all the frequency locations.
%         % ....
%         % ....
%         %%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
%         
%         rmax=find(data==Data_sort(j));
% %         if X(rmax)<th
% %             xmax=zeros(maxNumPeaks,1);
% %             ymax=zeros(maxNumPeaks,1);
% %             peaks_collect=[ymax xmax];
% %             break;
% %         else
%         if rmax==1
%             L_n= 0; L_m=data(rmax); L_p=data(rmax+1);
%             A=(L_p+L_n-2*(L_m))/2;
%             B=(L_p-L_n)/2;
%             C=L_m;
%             xmax=-B/2/A+F(rmax);%interpolate around the interested point
%             ymax=C-(B^2)/4/A;
%             
%         else
%             L_n= data(rmax-1); L_m=data(rmax); L_p=data(rmax+1);
%             A=(L_p+L_n-2*(L_m))/2;
%             B=(L_p-L_n)/2;
%             C=L_m;
%             xmax=-B/2/A+F(rmax);%interpolate around the interested point
%             ymax=C-(B^2)/4/A;
%             
%         end
%         peaks_collect=[peaks_collect;[ymax xmax]];
%         %%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
%         
%         
%     end
end
%% Return the list of amps and freqs in the order of ascending frequency
%%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
%     peaks_sort=sort(peaks_collect,1);
%     amps=peaks_sort(:,1);
%     freqs=peaks_sort(:,2);
%%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
        

end