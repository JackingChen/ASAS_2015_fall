ff=1:1:length(spectra1);
fs=500;
[psor,lsor] = findpeaks(spectra1);
hold off;;
plot(ff,spectra1);
peaks=zeros(length(spectra1),1);

count=1;
for i=1:length(spectra1)
    if(spectra(i)==lsor(count))
        print(lsor(count),i)
        peaks(i)=psor(count);
        count=count+1;
    end
end

hold on
plot(ff,peaks)