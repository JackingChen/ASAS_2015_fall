function out = ACF(wav_file,fs,window_length)

hop_size = window_length/2;
wav_length = length(wav_file);
st_num = ceil(wav_length/hop_size);
ed_num = st_num;
frame_num = st_num;
ACF_func = cell(frame_num,1);
st = zeros(st_num,1);
ed = zeros(ed_num,1);
out.pitch = zeros(frame_num,1);

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
        %%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%
        ACF_func{ii} = autocorr(wav_file(st(ii):ed(ii)),70);
%         xx=1:1:length(ACF_func{ii});
%         
%         if ii>5 && ii<10
%             figure()
%             plot(xx,ACF_func{ii})
%         end
        [pks loc]=findpeaks(ACF_func{ii},'MinPeakHeight',0.1);
%         pks_sarray=sort(pks,'descend');
%         largepks=pks(3);
%         T=loc(find(pks==largepks));
        T=loc(2);
        [V,I] = sort(ACF_func{ii},'descend');
        period=T/fs;
        pitch=1/period;
        out.pitch(ii) = pitch;
        %%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%
    else
        ACF_func{ii} = zeros(window_length,1);
        out.pitch(ii) = 0;
    end
end

out.ACF_func = ACF_func;

end