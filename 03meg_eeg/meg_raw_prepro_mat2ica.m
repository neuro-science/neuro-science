function errCode = meg_raw_prepro_mat2ica(p, f, is)
	% % % updated 20/09/2022 by wp for ica of clean data
	% % % copied from f04_tp_rawdata_ica
	% % % The first datasets were for TP4
	%% 1. initialize
	% % % 1.1 prepare outputs	
	errCode = 0;
	ps = [];

	% % % 1.2 check inputs	
	if nargin < 3
		errCode = sprintf( ...
			'Three parameters are needed: para, input and output');
		return;
	elseif ~exist(f.aus{is}, 'file')
		errCode = sprintf( ...
			'Input file %s was not found, exiting... \n', f.ein{is});
		return;
% 	elseif exist(f.skp{is}, 'file')
% 		errCode = sprintf( ...
% 			'The data %s is to skip, exiting... \n', f.idx{is});
% 		return;
	elseif exist(f.on{is}, 'file')
		errCode = sprintf( ...
			'The data %s is being processed, exiting... \n', f.idx{is});
		return;
	elseif exist(f.off{is}, 'file')
		errCode = sprintf( ...
			'The data %s was processed before, exiting... \n', f.idx{is});
		return;
% 	elseif exist(f.ica{is}, 'file')
% 		errCode = sprintf( ...
% 			'Ouput file %s exists, done before? exiting... \n', f.ica{is});
% 		return;
	end
	
	%% 2. prepare data files etc.
	% % % 2.1 place holder and welcome	
	save(f.on{is}, 'is');
	plt_myPrintLine([]);
	ttt = clock;
	fprintf('Process started on %s @%02d:%02d:%02.0f ...\n', ...
		f.idx{is}, ttt(4:6));
	c.fname = f.ica{is};

	% % % 2.2 load data 
	v = load(f.aus{is}, 'meg', 'eeg', 'hmv', 'ps');
	mdata = v.meg;
	mdata(v.ps.bad, :) = [];
	edata = v.eeg;
	edata(v.ps.bad, :) = [];
	ldata = v.hmv;
	ldata(v.ps.bad, :) = [];

	% % % 2.3 some info update	
	c.nPoints = size(mdata, 1);
	c.iTag = f.idx{is};
	nChansInfo = v.ps.nChans;
	c.nChs = size(mdata, 2);
	if abs(nChansInfo - c.nChs) > 0.1
		errCode = sprintf( ...
			'The number of channels mismatch in data and para!');
		return;
	end		
	ttt = clock;
	fprintf('Data loaded @%02d:%02d:%02.0f!\n', ttt(4:6));
	
	% % % 2.4 preapare locations	
	sChanID = 1 : c.nChs;
	senLocs = [v.ps.hdr.grad.coilpos(sChanID, :), ...
		v.ps.hdr.grad.coilori(sChanID, :)];
	loc2d = sensor3d2sensor2d(senLocs(:, 1:3));
	c.loc2d(:, 1) = -loc2d(:, 2);
	c.loc2d(:, 2) = loc2d(:, 1);	
	clear v;

	% % % 2.5 pick the ica waves position
	idx = linspace(0.1, 0.9, p.pre.icaWavePerData)';	%evenly spaced
	rand_shift = min(0.1, idx(2) - idx(1)) * (rand(1000, 1) - 0.5);
	idx = round((idx + rand_shift(1 : p.pre.icaWavePerData)) * c.nPoints);


	%% 3. do computation
	% % % 3.1 do PCA / select components	
	[D, E, infoRatio, nICs] = meg_raw_preFASTICA(mdata, ...
		p.pre.infoRatioICs, p.pre.numLim4ICs, p.pre.stepRatioICs);
	ttt = clock;
	fprintf('PCA computed @%02d:%02d:%02.0f!\n', ttt(4:6));

	% % % 3.2 do ICA
	[ica, A, W] = fastica(mdata', 'verbose', 'off', 'pcaE', E, 'pcaD', D);
	clear mdata;
	c.nICs = size(A, 2);
	ttt = clock;
	fprintf('ICA computed @%02d:%02d:%02.0f!\n', ttt(4:6));

	% % % 3.3. prepare saving the results for ICA
	c.A = A;
	c.W = W;
	c.D = D;
	c.E = E;
	c.infoRatio = infoRatio;

	% % % 4. prepare plot data
	% % % 4.1 topo data of ICA
	c.Z = zeros([c.nChs, c.nICs]);
	for ic = 1 : c.nICs
		c.Z(:, ic) = double(mean(A(:, ic) * ica(ic, :), 2));
	end

	% % % 4.2 prepare EOG/EKG data
	eog(:, 1) = edata(:, 2) - edata(:, 1);	% veog
	eog(:, 2) = edata(:, 4) - edata(:, 3);	% heog
	eog(:, 3) = edata(:, 5);					% ekg
	eog(:, 4) = mean(edata(:, 1:4), 2);		% reog
	clear edata;
	ica = ica';
	% % % correlation of EOG/EKG data
	[c.acc, c.acs] = meg_raw_corr_ICA_vs_EEG(ica, eog);
	% % % take the logorithm of p-value, 0.05-1.3; 0.01-2;
	c.acs = -log10(c.acs);										
	clear eog;

	% % % 4.3 ica data for curves
	c.xts = p.pre.icaWaveDispRange(1) : p.pre.icaWaveDispResolution : ...
		p.pre.icaWaveDispRange(2);
	c.nxts = numel(c.xts);
	c.yts = zeros(c.nxts, c.nICs, p.pre.icaWavePerData);
	c.pntWavePlot = bsxfun(@plus, c.xts', idx');
	for k = 1 : p.pre.icaWavePerData
		c.yts(:, :, k) = ica(c.pntWavePlot(:, k), :);
	end
	c.yts = permute(c.yts, [1 3 2]); %[pts, trl, ic]
	c.yts = bsxfun(@rdivide, c.yts, max(c.yts));

	% % % 4.4. plot of head motion ICA
	c.nm = floor(size(ica, 1) / p.pre.rawSampleRate / 60);
	ldata = bsxfun(@minus, ldata, nanmean(ldata));
	l2 = reshape(ldata(1 : p.pre.rawSampleRate * 60 * c.nm, :), ...
		p.pre.rawSampleRate * 60, c.nm, size(ldata, 2));
	td = squeeze(nanmean(l2, 1));
	c.trl = bsxfun(@rdivide, td, 20);	% 20 mm as base
	% % % correlation of head motion ICA
	[c.lcc, c.lcs] = meg_raw_corr_ICA_vs_EEG(ica, ldata);	% correlation coefficients/significants
	clear ldata td l2;

	% % % 4.5 spectra data for ICA
	tf0 = abs(fft(ica));	%time-frequency [npts, c.nICs]
	nn = ceil(c.nPoints / 2);	%number of frequencies
	ff = linspace(0, p.pre.rawSampleRate/2, nn + 1);	%frequencies
	ff(1) = [];
	% % % select range
	sid = find(ff <= p.pre.freqLowpassMEG & ff >= p.pre.freqHighpassMEG);	
	sid2 = round(linspace(sid(1), sid(end), p.pre.icaFreqDispNumber)); 
	yy = tf0(sid2, :);	% power for frequency
	MM = nanmax(yy(:));	% normalize base
	c.yfs = yy ./ MM;		% normalize to 0-1 
	c.xfs = ff(sid2);		% the xx
	clear tf0 nn sid MM ff yy MM sid2;

	% % % 4.6 trial evolution
	ica = reshape(ica(1 : p.pre.rawSampleRate * 60 * c.nm, :), p.pre.rawSampleRate * 60, c.nm, c.nICs);
	tf0 = abs(fft(ica));	%time-frequency [npts, c.nICs]
	ff = linspace(0, p.pre.rawSampleRate/2, p.pre.rawSampleRate * 30 + 1);	%frequencies
	ff(1) = [];
	fm2 = 20;
	sid = ff <= p.pre.freqLowpassMEG & ff >= fm2;	% select range
	td = squeeze(nanmean(tf0(sid, :, :), 1));
	c.trv = bsxfun(@rdivide, td, max(td));
	clear ff sid tf0 td fm2;

	%% 5. save the data
	save(f.ica{is}, 'c');
	ttt = clock;
	fprintf('ICA saved @%02d:%02d:%02.0f!\n', ttt(4:6));
	plt_myPrintLine([]);
	save(f.off{is}, 'is');

end

