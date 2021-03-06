function out = spectrum_analysis(wav_file,fs,window_length,poly_n)

hop_size = window_length/2;
wav_length = length(wav_file);
st_num = ceil(wav_length/hop_size);
ed_num = st_num;
frame_num = st_num;
spectra = cell(frame_num,1);
st = zeros(st_num,1);
ed = zeros(ed_num,1);
window_vec = hamming(window_length);

for ii = 1:frame_num
    if ii == 1
        st(ii) = 1 + (window_length)*(ii-1);
        ed(ii) = window_length*ii;
    elseif ii == frame_num
        st(ii) = st(ii-1) + hop_size;
        ed(ii) = wav_length;
    else
        st(ii) = st(ii-1) + hop_size;
        ed(ii) = ed(ii-1) + hop_size;
    end
end




for ii = 1:frame_num
    if ii ~= frame_num                
        temp = abs(fft(wav_file(st(ii):ed(ii)).*window_vec));
        spectra{ii} = zeros(floor(floor(length(temp)))/2+1,2);
        spectra{ii}(:,1) = 0:fs/length(temp):fs/2;
        spectra{ii}(:,2) = temp(1:floor(length(temp))/2+1);
    else
        spectra{ii} = zeros(floor(length(temp))/2+1,2);
        spectra{ii}(:,1) = 0:(fs/length(temp)):(fs/2);
    end
end

out.spectra = spectra;
out.pitch = zeros(frame_num,1);

starting=5;
ending=10;
count=0;

for ii = 1:frame_num-1
    %restrict the pitch 
    spectra{ii}(find(spectra{ii}(:,1)<50 | spectra{ii}(:,1)>500),2)=0;
    if poly_n>0
        % model the spectrum with ploynomial cn*x^poly_n+...c1*x^1+c0
        p=polyfit(spectra{ii}(:,1),spectra{ii}(:,2),poly_n); %train parameter p
        spectra{ii}(:,2)=polyval(p,spectra{ii}(:,1));% use p to reconstruct spectra
    end
    [pks,locs] = findpeaks(spectra{ii}(:,2));
    xx=1:1:length(spectra{ii}(:,2));
    count=ii;
    if (count<ending && count>starting) || (count<40 && count>35)
       figure
        plot(spectra{ii}(:,1),spectra{ii}(:,2),'o');
        title(num2str(ii))
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(pks)
        out.pitch(ii)=0;
    else
        [~,I] = max(pks);
        out.pitch(ii) = spectra{ii}(locs(I),1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%%%
end
end

