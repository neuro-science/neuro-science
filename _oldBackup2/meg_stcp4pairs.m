function [out1, out2] = meg_stcp4pairs(tf_data, pairs, mode, trial_id, time_idx, freq_idx, taperList, no_rebose_flag, no_stcp_flag)

% % % written 27/04/18 by wp for better documentation
% % % The variable tf_data shall be a cell array format {Time, Freq}
% % % The cell element shall be [channel, trial, taper]
% % % The varible pairs shale be 2-collumn array (ch1, ch2) 
% % % mode is a string which can be 'coh', 'plv' or 'imCoh'.
% % % stcp is a format of array [connection, trial]

	%% preparison
	tic;
	% % % check inputs and get info
	if nargin < 3
		error('At least 3 inputs needed:\n tf_data, pairs and mode!');
	else
		[nTs, nFs] = size(tf_data);
		[nChs, nTrls, tmp2] = size(tf_data{1, 1});
	end
	% % % optional paras
	if nargin < 4 || isempty(trial_id)
		trial_id = 1 : nTrls;
		nTrials = nTrls;
	else
		nTrials = length(trial_id);
	end
	
	if nargin < 5 || isempty(time_idx)
		time_idx = 1 : nTs;
	end
	numT = length(time_idx);

	if nargin < 6 || isempty(freq_idx)
		freq_idx = 1 : nFs;
	end
	numF = length(freq_idx);
	
	if nargin < 7 || isempty(taperList)
		nTapers = 1;
	else
		nTapers = reshape(taperList(freq_idx), [1 1 1 numF]);
	end
	
	if nargin < 8 || isempty(no_rebose_flag)
		no_rebose_flag = false;
	end
		
	if nargin < 9 || isempty(no_stcp_flag)
		no_stcp_flag = false;
	end
		
	% % % edge data info
	[nCons, tmp1] = size(pairs); %here cons mean connection
	maxCh = max(pairs(:));
	maxF = max(freq_idx);
	maxT = max(time_idx);
	maxTrials = max(trial_id);
	
	% % % check data
	if tmp1 ~= 2
		error('input pairs format error: not 4 collumns');
	end
	if maxCh > nChs
		error('number of channels error!');
	end
	if maxT > nTs
		error('number of time points error!');
	end
	if maxF > nFs
		error('number of frequencies error!');
	end
	if maxTrials > nTrls
		error('number of trials error!');
	end
	
	% % % innitiate connections
	pw = zeros(nCons, nTrials, numT, numF, 2);
	csd = zeros(nCons, nTrials, numT, numF);
	
	%% do data collection
	for it = 1 : numT
		for iq = 1 : numF
			d1 = tf_data{time_idx(it), freq_idx(iq)}(pairs(:, 1), trial_id, :);
			d2 = tf_data{time_idx(it), freq_idx(iq)}(pairs(:, 2), trial_id, :);
			pw(:, :, it, iq, 1) = mean(d1 .* conj(d1), 3);
			pw(:, :, it, iq, 2) = mean(d2 .* conj(d2), 3);
			csd(:, :, it, iq) = mean(d1 .* conj(d2), 3); %[con, trl]
			if ~no_rebose_flag
				fprintf('%04d-%02d-%02d %02d:%02d:%02d, done for frequency %d of %d @time %d of %d.\n', round(clock), iq, numF, it, numT);
			end
		end
	end
	fprintf('%04d-%02d-%02d %02d:%02d:%02d, data prepared.\n', round(clock));

	
	%% compute power
	ers = sum(pw, 5); 

	%% compute coherence
	switch lower(mode)
		case {'plv', 'p'}
			plv0 = csd ./ sqrt(pw(:, :, :, :, 1) .* pw(:, :, :, :, 2)); %[con, trl, t, f]
			s0 = sum(plv0, 2); %[con, 1, t, f]
			plva = abs(s0 / nTrials); %[con, 1, t, f]
			fprintf('%04d-%02d-%02d %02d:%02d:%02d, plv group done.\n', round(clock));
			plza = 23 / 20 * (sqrt(bsxfun(@times, log(1 - plva.^2), 2 - nTrials*nTapers*2)) - 23 / 20); %[con, 1, t, f]
			fprintf('%04d-%02d-%02d %02d:%02d:%02d, plv z-transformed.\n', round(clock));
			if no_stcp_flag
				out1 = plva;
				out2 = plza;
			else
				plvi = abs(bsxfun(@minus, s0, plv0) / (nTrials - 1)); %[con, trl, t, f]
				fprintf('%04d-%02d-%02d %02d:%02d:%02d, plv individual done.\n', round(clock));
				plzi = 23 / 20 * (sqrt(bsxfun(@times, log(1 - plvi.^2), 2 - (nTrials-1)*nTapers*2)) - 23 / 20); %[con, trl, t, f]
				fprintf('%04d-%02d-%02d %02d:%02d:%02d, plv individual z-transformed.\n', round(clock));
				stcp = bsxfun(@minus, nTrials * plza, (nTrials - 1) * plzi);
				fprintf('%04d-%02d-%02d %02d:%02d:%02d, stcp done.\n', round(clock));
				out1 = stcp;
				out2 = ers;
			end
			clear plv0 plva plvi plza plzi s0;
		case {'coh', 'c'}
			spw = sum(pw, 2);	%[con, 1, t, f, 2]
			coha = abs(sum(csd, 2)) ./ sqrt(spw(:, :, :, :, 1) .* spw(:, :, :, :, 2)); %[con, 1, t, f]
			coza = 23 / 20 * (sqrt(bsxfun(@times, log(1 - coha.^2), 2 - nTrials*nTapers*2)) - 23 / 20); %[con, 1, t, f]
			if no_stcp_flag
				out1 = coha;
				out2 = coza;
			else
				cohi = zeros(nCons, nTrials, numT, numF);
				for tr = 1 : nTrials
					spw1 = sum(pw(:, [1:tr-1, tr+1:end], :, :, :), 2);
					cohi(:, tr, :, :) = abs(sum(csd(:, [1:tr-1, tr+1:end], :, :), 2)) ./ sqrt(spw1(:, :, :, :, 1) .* spw1(:, :, :, :, 2)); %[con, 1, t, f]
				end
				cozi = 23 / 20 * (sqrt(bsxfun(@times,log(1 - cohi.^2), 2 - (nTrials-1)*nTapers*2)) - 23 / 20); %[con, trl, t, f]
				stcp = bsxfun(@minus, nTrials * coza, (nTrials - 1) * cozi);
				out1 = stcp;
				out2 = ers;
			end
			clear spw spw1 coha coza cohi cozi;
		case {'imcoh', 'i', 'imc'}
			spw = sum(pw, 2);	%[con, 1, t, f, 2]
			coha = imag(sum(csd, 2)) ./ sqrt(spw(:, :, :, :, 1) .* spw(:, :, :, :, 2)); %[con, 1, t, f]
			coza = 23 / 20 * (sqrt(bsxfun(@times, log(1 - coha.^2), 2 - nTrials*nTapers*2)) - 23 / 20); %[con, 1, t, f]
			if no_stcp_flag
				out1 = coha;
				out2 = coza;
			else
				cohi = zeros(nCons, nTrials, numT, numF);
				for tr = 1 : nTrials
					spw1 = sum(pw(:, [1:tr-1, tr+1:end], :, :, :), 2);
					cohi(:, tr) = imag(sum(csd(:, [1:tr-1, tr+1:end], :, :), 2)) ./ sqrt(spw1(:, :, :, :, 1) .* spw1(:, :, :, :, 2)); %[con, 1, t, f]
				end
				cozi = 23 / 20 * (sqrt(bsxfun(@times,log(1 - cohi.^2), 2 - (nTrials-1)*nTapers*2)) - 23 / 20); %[con, trl, t, f]
				stcp = bsxfun(@minus, nTrials * coza, (nTrials - 1) * cozi);
				out1 = stcp;
				out2 = ers;
			end
			clear spw spw1 coha coza cohi cozi;
		otherwise
			error('unkown coherence manner!');
	end
end %end of function