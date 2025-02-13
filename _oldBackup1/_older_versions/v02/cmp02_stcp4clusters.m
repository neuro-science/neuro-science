function [out1, out2] = cmp02_stcp4clusters(tf_data, edges, mode, trial_id, diff_time_idx, taperList, no_stcp_flag)

% % % updated 29/03/17 by wp, renewed with better comments
% % % updated 19/04/16 by wp, change from cmp02_Coupling4SelectPairs to current name
% % % updated 07/04/16 by wp for corrrection for a small error
% % % updated 26/06/14 by wp for better documentation

% % % stcp = wpa_0GEN_F101STCP(tf_data, edges)
% % % The variable tf_data shall be a cell array format {Time, Freq}
% % % The cell element shall be [channel, trial, taper]
% % % The varible edges shale be 4-collumn array (ch1, ch2, freq, time) 
% % % mode is a string which can be 'coh', 'plv' or 'imCoh'.
% % % stcp is a format of array [connection, trial]

	%% preparison
	tic;
	% % % check inputs and get info
	if nargin < 3
		error('At least 3 inputs needed:\n tf_data, edges and mode!');
	else
		[nTs, nFs] = size(tf_data);
		[nChs, nTrls, tmp2] = size(tf_data{1, 1});
	end
	% % % optional paras
	if nargin < 7 || isempty(no_stcp_flag)
		no_stcp_flag = false;
	end
	if nargin < 6 || isempty(taperList)
		nTapers = 1;
	else
		nTapers = taperList(edges(:, 3));
	end
	if nargin < 5 || isempty(diff_time_idx)
		diff_time_idx = 0;
	end
	if nargin < 4 || isempty(trial_id)
		trial_id = 1 : nTrls;
		nTrials = nTrls;
	else
		nTrials = length(trial_id);
	end
		
	% % % edge data info
	[nCons, tmp1] = size(edges); %here cons mean connection
	maxCh = max(max(edges(:, 1:2)));
	maxF = max(edges(:, 3));
	maxT = max(edges(:, 4));
	
	% % % check data
	if tmp1 ~= 4
		error('input edges format error: not 4 collumns');
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
	if nTrials > nTrls
		error('number of trials error!');
	end
	
	% % % innitiate connections
	pw = zeros(nCons, nTrials, 2);
	csd = zeros(nCons, nTrials);
	
	%% do data collection
	for c = 1 : nCons
		t = edges(c, 4) + diff_time_idx;
		f = edges(c, 3);
		ch1 = edges(c, 1);	%be careful on the sequence
		ch2 = edges(c, 2);
		
		d1 = tf_data{t, f}(ch1, trial_id, :);
		d2 = tf_data{t, f}(ch2, trial_id, :);
		csd(c, :) = mean(d1 .* conj(d2), 3); %[con, trl]
		pw(c, :, 1) = mean(d1 .* conj(d1), 3);	%power [con_node1, tr]
		pw(c, :, 2) = mean(d2 .* conj(d2), 3);	%power [con_node2, tr]
	end
	
	%% compute power
	ers = sum(pw, 3); 

	%% compute coherence
	switch lower(mode)
		case {'plv', 'p'}
			plv0 = csd ./ sqrt(pw(:, :, 1) .* pw(:, :, 2)); %[con, trl]
			s0 = sum(plv0, 2); %[con, 1]
			plva = abs(s0 / nTrials); %[con, 1]
			plza = 23 / 20 * (sqrt(log(1 - plva.^2) .* (2 - nTrials*nTapers*2)) - 23 / 20); %[con, 1]
			if no_stcp_flag
				out1 = plva;
				out2 = plza;
			else
				plvi = abs(bsxfun(@minus, s0, plv0) / (nTrials - 1)); %[con, trl]
				plzi = 23 / 20 * (sqrt(bsxfun(@times, log(1 - plvi.^2), 2 - (nTrials-1)*nTapers*2)) - 23 / 20); %[con, trl]
				stcp = bsxfun(@minus, nTrials * plza, (nTrials - 1) * plzi);
				out1 = stcp;
				out2 = ers;
			end
			clear plv0 plva plvi plza plzi s0;
		case {'coh', 'c'}
			spw = sum(pw, 2);	%[con, 1, 2]
			coha = abs(sum(csd, 2)) ./ sqrt(spw(:, :, 1) .* spw(:, :, 2)); %[con, 1]
			coza = 23 / 20 * (sqrt(log(1 - coha.^2) .* (2 - nTrials*nTapers*2)) - 23 / 20); %[con, 1]
			if no_stcp_flag
				out1 = coha;
				out2 = coza;
			else
				cohi = zeros(nCons, nTrials);
				parfor tr = 1 : nTrials
					spw1 = sum(pw(:, [1:tr-1, tr+1:end], :), 2);
					cohi(:, tr) = abs(sum(csd(:, [1:tr-1, tr+1:end]), 2)) ./ sqrt(spw1(:, :, 1) .* spw1(:, :, 2)); %[con, 1]
				end
				cozi = 23 / 20 * (sqrt(bsxfun(@times,log(1 - cohi.^2), 2 - (nTrials-1)*nTapers*2)) - 23 / 20); %[con, trl]
				stcp = bsxfun(@minus, nTrials * coza, (nTrials - 1) * cozi);
				out1 = stcp;
				out2 = ers;
			end
			clear spw spw1 coha coza cohi cozi;
		case {'imcoh', 'i', 'imc'}
			spw = sum(pw, 2);	%[con, 1, 2]
			coha = imag(sum(csd, 2)) ./ sqrt(spw(:, :, 1) .* spw(:, :, 2)); %[con, 1]
			coza = 23 / 20 * (sqrt(log(1 - coha.^2) .* (2 - nTrials*nTapers*2)) - 23 / 20); %[con, 1]
			if no_stcp_flag
				out1 = coha;
				out2 = coza;
			else
				cohi = zeros(nCons, nTrials);
				parfor tr = 1 : nTrials
					spw1 = sum(pw(:, [1:tr-1, tr+1:end], :), 2);
					cohi(:, tr) = imag(sum(csd(:, [1:tr-1, tr+1:end]), 2)) ./ sqrt(spw1(:, :, 1) .* spw1(:, :, 2)); %[con, 1]
				end
				cozi = 23 / 20 * (sqrt(bsxfun(@times,log(1 - cohi.^2), 2 - (nTrials-1)*nTapers*2)) - 23 / 20); %[con, trl]
				stcp = bsxfun(@minus, nTrials * coza, (nTrials - 1) * cozi);
				out1 = stcp;
				out2 = ers;
			end
			clear spw spw1 coha coza cohi cozi;
		otherwise
			error('unkown coherence manner!');
	end
end %end of function