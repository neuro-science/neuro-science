% % % To compute itpc needed phase information from raw data(timepoint x trial x channel)
% % % sRate in Hz, timeStart in ms, yFreqs in Hz, xTimes in ms, maxLength in ms.
% % % output preITPC (freq, time, chan, trial) 
function preITPC = meg_wavelet_preITPC (data, sRate, timeStart, xTimes, yFreqs, maxLength, maxNumCycle)

    %% 1. prepare the environment
    % 1.1 Check input parameters
    if nargin < 5
        error('You need to input data (timepoint x trial x channel), sampling rate, frequencies, times at least!');
    elseif nargin < 6 || isempty(maxLength)
        maxLength = 750; % Minimum 750 ms for window size as default
    elseif nargin < 7 || isempty(maxNumCycle)
        maxNumCycle = 6;
    end

    if min(yFreqs) < 1000 / maxLength
        fprintf('The max window size %d ms is not adequate for %d Hz', maxLength, min(yFreqs));
        return;
    end

    % 1.2 Convert time points from ms to sample indices
    xTimesInSamples = round((xTimes - timeStart) / 1000 * sRate) + 1;
    
    % 1.3 Preallocate preITPC array
    numFreqs = length(yFreqs);
    numTimes = length(xTimes);
    numTrials = size(data, 2);
    numChannels = size(data, 3);
    
    preITPC = zeros(numFreqs, numTimes, numChannels, numTrials);
    
    %% 2. do computation
    for fIdx = 1:numFreqs
        freq = yFreqs(fIdx);

        % Adjust the number of cycles based on frequency to ensure at least one cycle within the time window
        cycles = min(maxNumCycle, maxLength * freq / 1000); % Ensure at least one cycle
        
        % Generate wavelet
        s = cycles / (2 * pi * freq);
        t = -3*s : 1/sRate : 3*s;
        wavelet = exp(2 * 1i * pi * freq .* t) .* exp(-t.^2 / (2 * s^2));
        
        for chIdx = 1:numChannels
            for trialIdx = 1:numTrials
                % Convolve data with wavelet
                convolutionResult = conv(data(:, trialIdx, chIdx), wavelet, 'same');
                
                % Extract phase information at specified time points
                preITPC(fIdx, :, chIdx, trialIdx) = angle(convolutionResult(xTimesInSamples));
					 
					 clear convolutionResult;
            end
		  end
		  clear s t wavelet;
        fprintf('.');
    end
end
