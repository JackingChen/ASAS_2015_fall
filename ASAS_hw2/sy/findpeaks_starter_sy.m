function [amps,freqs]=findpeaks_starter_sy(X,freqq, maxNumPeaks)
% X=spectra{11}(:,2)
% freqq=spectra{11}(:,1)


% X is frequency component of one frame 
% maxNumPeaks is the number of peaks to be extracted

N=size(X(:),1);
data = 20*log10(abs(X));

data=data(:);
ind=find( (data>[data(1)-100;data((1:N-1)')]) ...
    & (data>=[data((2:N)'); data(N)-100]) );   % find peaks

Data=data(ind);
Data_sort=sortrows(Data);
Data_sort = Data_sort(end:-1:1);
Frame_PeakMaxSize = size(Data_sort,1); 

%% 狦 maxNumPeaks ⊿Τ禬筁–fame程peak计 :
if Frame_PeakMaxSize >= maxNumPeaks
    b=Data_sort(1:maxNumPeaks);
    peak=zeros(maxNumPeaks,1);
    freq=zeros(maxNumPeaks,1);
    for j=1:maxNumPeaks
        v=find(data==b(j));
        peak(j)=X(v);
        freq(j)=freqq(v);
    end
    
%% はぇ碞璶ㄏノず础猭ㄓ干peak
else    
    b=Data_sort(1:Frame_PeakMaxSize);
    Peak_lack = maxNumPeaks-Frame_PeakMaxSize; % Peak ぶ计
    aa=zeros(Frame_PeakMaxSize-1,1);
    
    for k=1:Frame_PeakMaxSize-1
        aa(k)=Data_sort(k+1)-Data_sort(k);
    end
    z=sort(aa,'descend');
    inter=zeros(Peak_lack,1);
    for bb=1:Peak_lack
        vv=find(aa==z(bb));
        inter(bb)=(Data_sort(vv)-Data_sort(vv+1))/2;
    end
    peak=zeros(maxNumPeaks,1);
    freq=zeros(maxNumPeaks,1);
    %%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%
    for j=1:Frame_PeakMaxSize
            %% Implement quadratic interpolation here, and save the
        %%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
        v=find(data==b(j));
        peak(j)=X(v);
        freq(j)=freqq(v);
        %%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%spectra{30}(:,2)
    end
    YI=interp1(data,freqq,inter);
    freq(Frame_PeakMaxSize+1:maxNumPeaks)=YI;
    peak(Frame_PeakMaxSize+1:maxNumPeaks)=10.^(inter./20);
end
%% Return the list of amps and freqs in the order of ascending frequency
    amps=peak;
    freqs=freq;

end
