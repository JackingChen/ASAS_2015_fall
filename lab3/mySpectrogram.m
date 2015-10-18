%   EE6641 ASAS Lab3: Spectrogram
%
%       mySpectrogram.m
%
%   Created by Hao-Wei Lee, 2015.09

function S = mySpectrogram( sig, window, Noverlap, Nfft )

    if nargin < 3, Noverlap = 0; end
    if nargin < 4, Nfft = length(window); end

    win_len = length(window);
    sig_len = length(sig);
    
    % Calculating # of frames ?N?:
    %   sig_len <= ( win_len - n_overlap ) * ( N - 1 ) + win_len
    %           <= ( win_len - n_overlap ) * N + n_overlap
    %   N       >= ( sig_len - n_overlap ) / ( win_len - n_overlap )
    Nframes = ceil( ( sig_len - Noverlap ) / ( win_len - Noverlap ) );
    
    % Padding the signal with zero s.t. ( sig_zp_len - n_overlap ) is
    % divisibly by ( win_len - n_overlap )
    sig_zp = zeros( (win_len - Noverlap) * Nframes + Noverlap, 1 );
    sig_zp( 1:sig_len ) = sig(:);
    
    % Initializing spectrogram matrix
    S = complex( zeros( Nfft, Nframes ) );
    
    % Generating signals frame-by-frame
    % Says:
    %   sig_len = 20
    %   win_len = 5
    %   n_overlap = 2
    %   
    %   start_idx = [ 1, 4, 7, 10, 13, 16 ]
    %   sig_idx_frames =
    %       [                                       [
    %           1, 4, 7, 10, 13, 16;                    0, 0, 0, 0, 0, 0;
    %           1, 4, 7, 10, 13, 16;                    1, 1, 1, 1, 1, 1;
    %           1, 4, 7, 10, 13, 16;        +           2, 2, 2, 2, 2, 2;
    %           1, 4, 7, 10, 13, 16;                    3, 3, 3, 3, 3, 3;
    %           1, 4, 7, 10, 13, 16;                    4, 4, 4, 4, 4, 4;
    %                                  ]                                    ]
    %
    %       [
    %           1, 4,  7, 10, 13, 16;
    %           2, 5,  8, 11, 14, 17;
    %   =       3, 6,  9, 12, 15, 18;
    %           4, 7, 10, 13, 16, 19;
    %           5, 8, 11, 14, 17, 20;
    %                                   ]
    %
    %   sig_frams =
    %       [
    %           sig_zp[1], sig_zp[4], sig_zp[7],  ...
    %           sig_zp[2], sig_zp[5], sig_zp[8],  ...
    %              .          .           .       ...
    %              .          .           .       ...
    %           sig_zp[5], sig_zp[8], sig_zp[11], ...
    %                                                   ]
    %
    start_idx = 1 : (win_len - Noverlap) : (length(sig_zp) - win_len + 1);%start index= 1: interval of each step: 
    sig_idx_frames = repmat( start_idx, [win_len 1] ) + repmat( (1:win_len)' - 1, [1 length(start_idx)] );
    sig_frames = sig_zp( sig_idx_frames );
    

    
    %%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%
    
    
    
    for col = 1:size(sig_frames, 2)
        % get the signal of the col'th frame
        %   dummy = ...
        dummy = sig_frames( :, col );
        
        % dot-multiply by the window
        %   dummy = ...
        dummy = dummy .* window;
        
        % compute the spectrum and store in the col's frame
        S( :, col ) = fft( dummy, Nfft );
    end
    
    
    
    %%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%
    
    S = S( 1:floor(end/2), : ); %take half of the spectrum
    
end

