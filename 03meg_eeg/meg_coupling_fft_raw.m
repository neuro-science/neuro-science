
function out = meg_coupling_fft_raw (data, pnts, trls, type)
% % % written 19/12/24, different trials in each cell
% % % `data` shall be points x trial x chans
% % % `pnts` is the logical index or range for the selected time span
% % % `trls` is the logical index or range for the selected trials
	
	% Parameters
	selected_data = data(pnts, trls, :); % Extract data for the selected time and trials
	[n_time, ~, ~] = size(selected_data);

	% Compute FFT across time for each trial and channel
	Y = fft(selected_data, [], 1); % FFT along the time dimension
	Y = Y(1:ceil(n_time/2), :, :); % Take only positive frequencies (Nyquist)
	
	% Compute auto-spectral densities for all channels
	Sxx = sum(abs(Y).^2, 2); % [freq x chan]

	% Compute cross-spectral densities for all channel pairs
	Cxy = pagemtimes(permute(Y, [3 2 1]), 'none', conj(permute(Y, [2 3 1])), 'none');

	% Compute coherency
	coherency = Cxy ./ sqrt(pagemtimes(permute(Sxx, [3 2 1]), permute(Sxx, [2 3 1]))); % [freq x chan x chan]

	% Extract the imaginary part of coherency
	if nargin > 3 && strncmp(type, 'i', 1)
		out = imag(coherency); % [freq x chan x chan]
	else
		out = abs(coherency); % [freq x chan x chan]
	end			

end %end of function