function [out1, out2] = meg_stcp4pac(src, ed, LF, HF, winsize, xtimes, xts, trl)

% % % updated 22/05/18 by wp, first draft based on meg_stcp4clusters

	%% preparison
	tic;
	% % % check inputs and get info
	if nargin < 7
		error('At least 4 inputs needed:\n 1src_data and 2ed, 3LF and 4HF, 5winsize, 6xtimes, 7xts!');
	elseif nargin < 8 || isempty(trl)
		trl = 1 : size(src, 3);
	end
	
% % % 	
% % % 		if nargin < 8 || isempty(st)
% % % 			st = 1 : npt;
% % % 		end
% % % 		
% % % 		end
% % % 		if nargin < 6 || isempty(ch2)
% % % 			ch2 = 1 : nv;
% % % 		end
% % % 		if nargin < 5 || isempty(ch1)
% % % 			ch1 = 1 : nv;
% % % 		end
% % % 		if isvector(st)
% % % 			d1 = reshape(src(ch1, st, trl), [numel(ch1), numel(st)*numel(trl)])';	%ch, point
% % % 			d2 = reshape(src(ch2, st, trl), [numel(ch2), numel(st)*numel(trl)])';
% % % 			clear src;
% % % 		else
% % % 			sz = size(st);
% % % 			if numel(sz) > 2
% % % 				error('Only one time series (2-d time segment) is allowed!');
% % % 			else
% % % 				%[ch, point, time bin, trials] -> [point, trials, ch, time bin]
% % % 				d1 = permute(reshape(src(ch1, st, trl), [numel(ch1), sz(1), sz(2), numel(trl)]), [2 4 1 3]);	
% % % 				d2 = permute(reshape(src(ch2, st, trl), [numel(ch2), sz(1), sz(2), numel(trl)]), [2 4 1 3]);	
% % % 			end
% % % 		end
% % % 	end
	
	%% do data collection
	for c = 1 : nCons
		t = ed(c, 4) + diff_time_idx;
		f = ed(c, 3);
		ch1 = ed(c, 1);	%be careful on the sequence
		ch2 = ed(c, 2);
		
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
	
	

%% x.15 pac of induced data 2018-05-22
function r = pac_seed2all_induced_con2(p, s, r, fpath)
					
	% % % check the existence of previous work
	r.pacName = [fpath, '3pros/05pac/', r.prefix, '.mat'];
	if exist(r.pacName, 'file')
		fprintf('We have this done before, skipped!\n');
		return;
	else
		save(r.pacName, 'r');
	end

	% % % read para data and reset
	% % % paras
	r.nbins = 18;
	r.nrepeats = 100;
	r.winsize = 250;
	r.fOrfer = 4;
	r.n1 = length(r.r1);
	r.n2 = length(r.r2);
	[r.fB1, r.fA1] = butter(r.fOrfer, 2 * [r.f1(1) r.f1(2)] / p.srate); 
	[r.fB2, r.fA2] = butter(r.fOrfer, 2 * [r.f2(1) r.f2(2)] / p.srate); 

	% % % do data in the loop for subjects
	for is = 1 : p.nsbjs
		if exist([fpath, '3pros/99tmp/', s{is}.datID, r.prefix, '_T', num2str(p.nts), '.mat'], 'file')
			fprintf('The subject was processed before, skipped!\n===\n');
			continue;
		end
		% % % prepare data for coupling
		fprintf('\n====\nTrying data %s:\n', s{is}.datID);
		fn_src = [fpath, s{is}.fn.src];
		s{is}.fn.tmp = [fpath, '3pros/99tmp/', s{is}.datID, r.prefix, '_T', num2str(p.nts), '.mat'];
		if exist(s{is}.fn.tmp, 'file')
			fprintf('Data were processed before, skipped!\n===\n');
			continue;
		elseif exist(fn_src, 'file')
			v = load(fn_src, 'src');
			sz = size(v.src);
			if (sz(1) - p.nv) || (sz(2) - p.npts)
				fprintf('Data size is odd, check it!');
				return;
			end
			src = v.src(:, :, s{is}.trls2(:));
			clear v;
			fprintf('Source data exist, erp removed...\n');
		else
			fprintf('Source data were missing, skipped!\n===\n');
			continue;
		end
		% % % work in loop for all times to get pac
		for it = 1 : p.nts
			% % % prepare file
			s{is}.fn.tmp = [fpath, '3pros/99tmp/', s{is}.datID, r.prefix, '_T', num2str(it), '.mat'];
			if exist(s{is}.fn.tmp, 'file')
				fprintf('Data were processed before, skipped!\n===\n');
				continue;
			else
				save(s{is}.fn.tmp, 'is');
				fprintf('Data will be processed...\n');
			end
			% % % para set
			r.tid = p.xtime <= p.xts(it) + r.winsize/2 & p.xtime >= p.xts(it) - r.winsize/2;
			r.m = length(find(r.tid));
			% % % prepare data for coupling
			src_now = reshape(src(:, r.tid, :), p.nv, r.m * s{is}.ntrls2*2);
			% % % LF angle
			d01 = src_now(r.r1, :)';
			tmp = [d01(1:1000, :); d01; d01(1:1000, :)];
			d11 = neu_preCrossFrequency(tmp, r.fA1, r.fB1, 'a');
			d11([1 : 1000, end-999 : end], :) = [];
			clear tmp;
			% % % HF Amplitude
			d02 = src_now(r.r2, :)';
			tmp = [d02(1:1000, :); d02; d02(1:1000, :)];
			d22 = neu_preCrossFrequency(tmp, r.fA2, r.fB2, 'A');
			d22([1 : 1000, end-999 : end], :) = [];
			clear tmp;
			% % % compute MI
			[pac1, pac_r1] = neu_crossfreq_MI(d11(1 : r.m * s{is}.ntrls2, :), d22(1 : r.m * s{is}.ntrls2, :), r.nbins, r.m, r.nrepeats);
			fprintf('Result1 computed on %04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
			[pac2, pac_r2] = neu_crossfreq_MI(d11(r.m * s{is}.ntrls2 + 1 : end, :), d22(r.m * s{is}.ntrls2 + 1 : end, :), r.nbins, r.m, r.nrepeats);
			fprintf('Result2 computed on %04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
			% % % save data
			save(s{is}.fn.tmp, 'pac1', 'pac_r1', 'pac2', 'pac_r2', 'r'); %%%!!!The data were deleted later	after summary
			fprintf('Results saved on %04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
			clear pac1 pac_r1 pac2 pac_r2;
		end
		clear src;
	end

	% % % summarize pac to one file - initialize
	pac = zeros(r.n1, r.n2, 2, p.nts, p.nsbjs, 'single');
	pac_r = zeros(r.n1, r.n2, r.nrepeats, 2, p.nts, p.nsbjs, 'single');
	% % % summarize pac to one file - loop on time and subjects
	for it = 1 : p.nts
		for is = 1 : p.nsbjs
			s{is}.fn.tmp = [fpath, '3pros/99tmp/', s{is}.datID,  r.prefix, '_T', num2str(it), '.mat'];
			v = load(s{is}.fn.tmp);
			pac(:, :, 1, it, is) = single(v.pac1);
			pac_r(:, :, :, 1, it, is) = single(v.pac_r1);
			pac(:, :, 2, it, is) = single(v.pac2);
			pac_r(:, :, :, 2, it, is) = single(v.pac_r2);
			clear v;
			fprintf('Data collected for subject %02d time %02d on %04d-%02d-%02d %02d:%02d:%02d.\n', is, it, round(clock));
		end
	end
	save(r.pacName, 'pac', 'pac_r', 'r');
	clear pac pac_r;
end
	
	
	
	
end %end of function