function [ output_args ] = plot_indicate( spectra )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ff=1:1:length(spectra);

[psor,lsor] = findpeaks(spectra);
hold off;;
plot(ff,spectra);
peaks=zeros(length(spectra),1);

count=1;
for i=1:length(spectra)
    if(count>length(psor))
            break
    end
    if(spectra(i)==psor(count))        
        peaks(i)=psor(count);
        count=count+1
    end
end

hold on
plot(ff,peaks,'o')

end

