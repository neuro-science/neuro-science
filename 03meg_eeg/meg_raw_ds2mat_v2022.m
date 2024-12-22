function errCode = meg_raw_prepro_ds2mat(p, f, is)
% % % written 18/09/2022 as summary of all preprocess before ICA
% % % The first datasets were for TP4
	%% 1. initialize
	% % % 1.1 prepare outputs	
	errCode = 0;
	ps = [];

	% % % 1.2 check inputs	
	if nargin < 3
		errCode = sprintf( ...
			'Three parameters are needed: para, fileInfo and ID');
		return;
	elseif ~exist(f.ein{is}, 'file')
		errCode = sprintf( ...
			'Input file %s was not found, exiting... \n', f.ein{is});
		return;
	elseif exist(f.on{is}, 'file')
		errCode = sprintf('The data is being processed, exiting... \n');
		return;
	elseif exist(f.off{is}, 'file')
		errCode = sprintf('The data was processed before, exiting... \n');
		return;
% 	elseif exist(f.aus{is}, 'file')
% 		errCode = sprintf( ...
% 			'Ouput file %s exists, done before? exiting... \n', f.aus{is});
% 		return;
	end
	
	%% 2. read in from raw data	
	% % % 2.1 place holder and welcome	
	save(f.on{is}, 'is');
	plt_myPrintLine([]);
	ttt = clock;
	fprintf('Process started on %s @%02d:%02d:%02.0f ...\n', ...
		f.idx{is}, ttt(4:6));

	% % % 2.2 read header and events	
	hdr = ft_read_header(f.ein{is});
	evt = ft_read_event(f.ein{is});

	% % % 2.3 read meg data	ps.nPoints
	megIDs = find(strncmp('M', hdr.label, 1));
	meg = ft_read_data(f.ein{is}, 'chanindx', megIDs); 
	sz = size(meg);
	meg = reshape(meg, sz(1), [])';
	[ps.nPoints, ps.nChans] = size(meg);

	% % % EEG read
	ps.nEEGChans = length(p.pre.eegLabels);
	[tmp0, tmp1, tmp2] = intersect(p.pre.eegLabels, hdr.label);
	[Y, I] = sort(tmp1);
	ps.eegIDs = tmp2(I);
	eeg = ft_read_data(f.ein{is}, 'chanindx', ps.eegIDs); 
	eeg = reshape(eeg, [ps.nEEGChans, ps.nPoints])';

	% % % HMV read
	ps.nHMVChans = length(p.pre.hmvLabels);
	[tmp0, tmp1, tmp2] = intersect(p.pre.hmvLabels, hdr.label);
	[Y, I] = sort(tmp1);
	ps.hmvIDs = tmp2(I);
	hmv = ft_read_data(f.ein{is}, 'chanindx', ps.hmvIDs); 
	hmv = reshape(hmv, [ps.nHMVChans, ps.nPoints])';
	clear Y I tmp*;
		
	% % % get paras	
	ps.evt = evt;
	ps.hdr = hdr;
	ttt = clock;
	fprintf('Loaded on %s @%02d:%02d:%02.0f ...\n', f.idx{is}, ttt(4:6));
	clear evt hdr;

	%% 3. clean up	
	% % % 3.1 remove breaks
	% % % get the segments
	ps.interBlockSegments(:, 1) = [1, f.endPoints{is} + ...
		round(p.pre.hinBlockInSec * p.pre.rawSampleRate)];
	ps.interBlockSegments(:, 2) = [f.startPoints{is} - ...
		round(p.pre.vorBlockInSec * p.pre.rawSampleRate), ps.nPoints];
	ps.interBlockSegments(ps.interBlockSegments < 1) = 1;
	ps.interBlockSegments(ps.interBlockSegments > ps.nPoints) = ps.nPoints;
	ps.interBlockSegments(ps.interBlockSegments(:, 2) <= ...
		ps.interBlockSegments(:, 1), :) = [];
	% % % remove the inter-block segments
	[meg, ps.ib] = meg_raw_rm_breaks(meg, ps.interBlockSegments);
	
	% % % 3.2 remove cliffs
	[meg, ps.zc] = meg_raw_rm_cliff(meg, p.pre.zThreshBaseCliff, ...
		p.pre.zThreshRatioCliff, p.pre.zRangeCliff);

	% % % 3.3 detrend to ease peak removing	
	meg1 = detrend(meg);

	% % % 3.4 remove peaks 1	
	X = 0; Y = 9e4; Z = p.pre.zThreshPeak; R = p.pre.zRangePeak;
	while Y > ps.nPoints / p.pre.zSegmentsNumThresh
		[meg, ps.zp1] = meg_raw_rm_peak(meg1 * p.pre.megScaleFactor, Z, R);
		X = X + 1;
		if isfield(ps.zp1, 'theSegments')
			Y = size(ps.zp1.theSegments, 1);
		else
			Y = 0;
		end
		fprintf('[%5d] segments @Z%5.2f @R%3d, counter %02d ...\n', ...
			Y, Z, R, X);
		Z = Z + p.pre.zThreshPeakStep;
		R = R + p.pre.zRangePeakStep;
	end
	clear meg1;
	ttt = clock;
	fprintf('Peak last removed @%02d:%02d:%02.0f!\n\n', ttt(4:6));
	
	% % % 3.5 check machine noise - to be improved
	% % % 		% % % 5.1.3 fft and save noisy scales // not 2nd of v2 again
	% % % 		meg = reshape(meg, [sz(2) sz(3) sz(1)]);	%[pt, tr, ch]
	% % % 		y = abs(fft(meg));
	% % % 		tfd = y(2:sz(2)/2+1, :, :);	%[f, tr, ch]
	% % % 		clear y;
	% % % 		dd = bsxfun(@eq, tfd, max(tfd)) & tfd > 0;
	% % % 		dd(sFreq, :, :) = 0;
	% % % 		r = permute(any(dd, 1), [2 3 1]);
	% % % 
	% % % 		ttt = clock;
	% % % 		if sum(r(:)) > sz(1) * sz(3) * 0.001
	% % % 			d = bsxfun(@rdivide, tfd, tfd(50, :, :));%normalized within trials to 50Hz
	% % % 			d3 = mean(d, 3, 'omitnan'); %average across channels
	% % % 			d2 = mean(d, 2, 'omitnan'); %average across trials
	% % % 			e3 = std(d, [], 3, 'omitnan'); %average across channels
	% % % 			e2 = std(d, [], 2, 'omitnan'); %average across trials
	% % % 			z = bsxfun(@rdivide,	bsxfun(@minus, d, d2), e2);
	% % % 			mz = mean(z, 3, 'omitnan');
	% % % 			ez = std(z, [], 3, 'omitnan');
	% % % 			save(f.tmp{is}, 'd', 'z', 'd2', 'e2', 'd3', 'e3', 'mz', 'ez', 'r');
	% % % 			fprintf('Dirty @%02d:%02d:%02.0f!\n\n', ttt(4:6));
	% % % 			clear d z d2 e2 d3 e3 mz ez dd tfd r;
	% % % 		else
	% % % 			save(f.tmp{is}, 'r');
	% % % 			fprintf('Clean @%02d:%02d:%02.0f!\n\n', ttt(4:6));
	% % % 			clear dd tfd r;
	% % % 		end
	% % % 		meg = reshape(meg, [sz(2)*sz(3) sz(1)]);	%[pt, tr, ch]

	% % % 3.6 filter meg
	sz = size(meg);
	if sz(1) < p.pre.fltSeparatePointThresh || ...
			sz(2) <= p.pre.fltSeparateChanThresh
		fprintf('Normal Filtering will be performed.\n')
		[meg, ps.zf] = meg_raw_filter(meg, p.pre.rawSampleRate, ...
			p.pre.freqHighpassMEG, p.pre.freqLowpassMEG, p.pre.fOrderMEG);
	else
		fprintf('Splitted Filtering will be performed.\n')
		meg1 = meg(:, 1 : p.pre.fltSeparateChanThresh);
		meg2 = meg(:, p.pre.fltSeparateChanThresh + 1 : end);
		[meg1, ps.zf] = meg_raw_filter(meg1, p.pre.rawSampleRate, ...
			p.pre.freqHighpassMEG, p.pre.freqLowpassMEG, p.pre.fOrderMEG);
		[meg2, ps.zf] = meg_raw_filter(meg2, p.pre.rawSampleRate, ...
			p.pre.freqHighpassMEG, p.pre.freqLowpassMEG, p.pre.fOrderMEG);
		meg = cat(2, meg1, meg2);
		clear meg1 meg2;
		if max(abs(size(meg) - sz)) > 0.001
			errCode = sprintf( 'XXX Data size mismatch! XXX \n');
			return;
		end
	end
	ttt = clock;
	fprintf('Filtered @%02d:%02d:%02.0f!\n\n', ttt(4:6));

	
	% % % 3.7 remove peaks 2	
	X = 0; Y = 9e4; Z = p.pre.zThreshPeak; R = p.pre.zRangePeak;
	while Y > ps.nPoints / p.pre.zSegmentsNumThresh
		[meg1, ps.zp2] = meg_raw_rm_peak(meg, Z, R);
		X = X + 1;
		if isfield(ps.zp2, 'theSegments')
			Y = size(ps.zp2.theSegments, 1);
		else
			Y = 0;
		end
		fprintf('[%5d] segments @Z%5.2f @R%3d, counter %02d ...\n', ...
			Y, Z, R, X);
		Z = Z + p.pre.zThreshPeakStep;
		R = R + p.pre.zRangePeakStep;
	end
	meg = meg1;
	clear meg1;
	ttt = clock;
	fprintf('Peak last removed @%02d:%02d:%02.0f!\n\n', ttt(4:6));

	% % % 3.8 filter eeg
	eeg = meg_raw_filter(eeg, p.pre.rawSampleRate, ...
		p.pre.freqHighpassMEG, p.pre.freqLowpassMEG, p.pre.fOrderMEG);
	ttt = clock;
	fprintf('EEG filtered @%02d:%02d:%02.0f!\n\n', ttt(4:6));

	% % % 3.9 summary of bad indices
	if isfield(ps.ib, 'segments')
		segs = ps.ib.segments;
	else
		segs = [];
	end
	if isfield(ps.zc, 'theSegments')
		segs = [segs; ps.zc.theSegments];
	end
	if isfield(ps.zp1, 'theSegments')
		segs = [segs; ps.zp1.theSegments];
	end
	if isfield(ps.zp2, 'theSegments')
		segs = [segs; ps.zp2.theSegments];
	end
	ps.bad = [];
	for k = 1 : size(segs, 1)
		ps.bad = [ps.bad, segs(k, 1) : segs(k, 2)];
	end
	ps.bad = unique(ps.bad);
	clear segs k;
	ttt = clock;
	fprintf('Preprocessed @%02d:%02d:%02.0f!\n\n', ttt(4:6));

	% % % 3.10 remove 50Hz meg/eeg
	[meg, ps.r50] = meg_raw_rm_50Hz(meg, p.pre.segLength50, ...
		ps.hdr.Fs, p.pre.lineNoiseFreq, ps.bad, ...
		[p.pre.freqHighpassMEG, p.pre.freqLowpassMEG]);
	[eeg, tmp] = meg_raw_rm_50Hz(eeg, p.pre.segLength50, ...
		ps.hdr.Fs, p.pre.lineNoiseFreq, ps.bad, ...
		[p.pre.freqHighpassMEG, p.pre.freqLowpassMEG]);

	% % % 3.11 clean up + save
	save(f.aus{is}, 'meg', 'eeg', 'hmv', 'ps', '-v7.3');
	clear ps meg eeg hmv sz;
	save(f.off{is}, 'is');
	ttt = clock;
	fprintf('Data saved @%02d:%02d:%02.0f!\n======\n\n', ttt(4:6));

end

