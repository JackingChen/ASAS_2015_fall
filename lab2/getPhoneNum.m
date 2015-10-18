function [ number ] = getPhoneNum( WAV_PATH )
%
%   EE6641 Lab2
%
%
%   EE6641 Lab2
%


Nnums = 10;
useFFT = 0;
filename =  WAV_PATH ;

%% Import the signal
[sig fs] = wavread(filename);
sig = sig(:, 1);
sig = sig(1 : floor(length(sig)/Nnums) * Nnums); %make sure it's multiple of 10(Nnums)

fprintf(['The signal length: ' num2str( length(sig)/fs ) ' sec.\nThe sampling frequency: ' num2str(fs) ' Hz.\n']);

digSigs = reshape( sig, [ length(sig)/Nnums, Nnums ] );
if length(sig) >4000
    spec_shorten=1
end
%make sure sig have same length
% digSigs=digSigs(1:floor(length(digSigs(:,1))/2), :);

tic;
if ~useFFT
    %% Construct DFT Matrix
    Ndft = size(digSigs, 1);
    %%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%
    para=-2*pi*j/Ndft;
    expo=exp(para);
    DFTMat = zeros( Ndft );
    for i =1:Ndft
        for l=1:Ndft
            DFTMat(i,l)=expo.^((i-1)*(l-1));
        end
    end
    %%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%
    

    
    %% Compute the spectrum of each piece of signal
    spectra = DFTMat * digSigs;
else
    Ndft = size(digSigs, 1);
    spectra = complex( zeros( Ndft, 10 ) );
    for ii = 1:10
        spectra(:, ii) = fft(digSigs(:, ii), Ndft);
    end
end
elapseTimeMillis = toc;

%% Post-processing the spectra
spectra = abs(spectra);
spectra = spectra(1:floor(length(spectra(:,1))/2), :);%take half of the wave
% my math, cos(?t+?)=12(e+i?t+e?i?t) there will be +? and -? be sensed,
% in addition, spectrum are replicated in order of mulitiple of fundamental
% frequency, so the frequency component is mirrored



%% Visualize the spectra
figure(1);
ff = ( 1:size(spectra, 1) )' / size(spectra, 1) * (fs/2);%divide the sample freq by the total points, each point represent certain part of freq

for ii = 1:10
    subplot( 10, 1, ii );
    plot( ff, spectra(:, ii) );
    xlim([0 2000]);
end

%% Figuring out what digit has been pressed
BTN_LIST = [
    '1', '2', '3', 'a'; ...
    '4', '5', '6', 'b'; ...
    '7', '8', '9', 'c'; ...
    '.', '0', '#', 'd' 
];
FREQ_LIST_VERT = [ 697, 770, 852, 941 ];
FREQ_LIST_HORIZ = [ 1209, 1336, 1477, 1633 ];



%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%
%get the freq
ii=1

for ii = 1:length(1:size(spectra, 2))

    spectra1=spectra(:,ii);
%     if spec_shorten==1
%     spectra1= spectra1(1:floor(length(spectra1(:,1))/2));
%     end
    %find two peaks
    %     ff=1:1:length(spectra1)
    %     figure(2)
    %     plot( ff, spectra1 );
    %     xlim([400 600]);
    threshold_peak=300;
    
    %     [pks locs]= findpeaks(spectra1);
    
    
    [pks loc]=findpeaks(spectra1);
    phase=0;
%     windowsize=5
%     numstep=floor(length(pks)/windowsize)
    for i =1:length(pks)
%         p_value=max(pks(1+(i-1)*windowsize:i*windowsize))
        if (pks(i)>threshold_peak && phase==0)
            freq_low=pks(i);
            phase=1;
        elseif (pks(i)>threshold_peak && phase==1)
            freq_high=pks(i);
            phase=0;
        else
            continue
        end
    end
    
    freq_pointl=loc(find(pks==freq_low));
    freq_pointh=loc(find(pks==freq_high));
    
    freql=ff(freq_pointl);
    freqh=ff(freq_pointh);
% figure(2)%******************    
% plot_indicate(spectra1)%****************** 
    
    
%sort the frequency
    if freql>freqh
        freq1(ii)=freqh;
        freq2(ii)=freql;
    elseif freqh>freql
        freq1(ii)=freql;
        freq2(ii)=freqh;
    else
        freq1(ii)=freql;
        freq2(ii)=freqh;
    end
end
    
% match the found peak to the number    
    threshold_HORIZ=50;
    threshold_VERT=30;
    
    
    
    %find number
    for i= 1:length(freq1)
        %check index on FREQ_LIST_VERT
        for index = 1:length(FREQ_LIST_VERT)
            if abs(FREQ_LIST_VERT(index)-freq1(i))<threshold_VERT
                VERT=index;
            else
                continue
            end
        end
        %check index on FREQ_LIST_HORTZ
        for index = 1:length(FREQ_LIST_HORIZ)
            if abs(FREQ_LIST_HORIZ(index)-freq2(i))<threshold_HORIZ
                HORIZ=index;
            else
                continue
            end
        end
        num(i)=BTN_LIST(VERT,HORIZ);
    end


number = num;

%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%



%% Print the result
if useFFT, method = 'fft';
else method = 'matrix dft'; end

fprintf( ['The phone number is: ' number '.\n'] );
fprintf( ['The elapsed time of Fourier Transform(' method '): ' num2str(elapseTimeMillis) ' ms.\n'] );

clear BTN_LIST FREQ_LIST_VERT FREQ_LIST_HORIZ ff filename ii method useFFT sig;

end
