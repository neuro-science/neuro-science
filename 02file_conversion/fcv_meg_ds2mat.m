function [p, hdr, evt, eeg, hmv, meg] = fcv_meg_ds2mat ...
	(inName, outName, flagWetRun, eegID, hmvID)
% % % Read the raw CTF MEG data and epoch them to trials
% % % 03/10/23	updated by wp: minor
% % % 27/03/23	re-written by wp: improve save for dry: no save empty filename
% % % 03/02/23	re-written by wp: add half-dry option
% % % 	flagWetRun: =0 [false] - read info without data
% % % 	flagWetRun: >0 - read data and save
% % % 	flagWetRun: <0 - read data but not save
% % % 08/03/22	re-written by wp: add dry option
% % % 17/07/14	re-written by wp: 
% % % fieldtrip componets <ft_read_header, ft_read_event, ft_read_data > were used
% % % converted from <wpa_0GEN_F01ReadCTFData>
% % % updated 04/09/2013 remove sorting part
% % % updated 03/09/2013 to make it more general
% % % updated 16/04/2013 to fix the bug of immediate start of experiment

	%% 1. check inputs
% 	ttt = clock;
% 	fprintf('Trying on data readin of %s @%02d:%02d:%02.0f...\n', inName, ttt(4:6));
	if nargin < 1 || isempty(inName)
		error('You have to specify the input file name!');
	elseif nargin < 2
		outName = [];
		eegID = [];
		hmvID = [];
		flagWetRun = false;
		fprintf('This is a dry read in without data!\n');
	elseif nargin < 3 || isempty(flagWetRun) || ~flagWetRun
		eegID = [];
		hmvID = [];
		flagWetRun = false;
		fprintf('This is a dry read in without data!\n');
	elseif nargin < 4
		eegID = [];
		hmvID = [];
		fprintf('no EEG/HM channels specified, only MEG data read!\n');
	elseif nargin < 5
		hmvID = [];
		fprintf('no HM channels specified, only MEG/EEG data read!\n');
	end
	
	%%	2. read header and event data
	hdr = ft_read_header(inName);	
	evt = ft_read_event(inName);
	
	%% 3. check info and get para set
	idx = strncmp('M', hdr.label, 1);
	p.megIDs = find(idx);
	p.megLabels = hdr.label(idx);
	p.numChannels = numel(p.megIDs);
	p.rawNumTrials = hdr.nTrials;
	p.rawPointsPerTrial = hdr.nSamples;
	p.rawSampleRate = hdr.Fs;
	p.rawNumPoints = p.rawPointsPerTrial * p.rawNumTrials;
	fs = {'p', 'hdr', 'evt'};
	
	%% 4. EEG read
	if ~isempty(eegID) && flagWetRun
		p.nEEGChans = length(eegID);
		if isnumeric(eegID)
			p.eegIDs = eegID;
			p.eegLabels	= hdr.label(eegID);
		elseif iscell(eegID)
			p.eegLabels	= eegID;
			[tmp0, tmp1, tmp2] = intersect(eegID, hdr.label);
			[Y, I] = sort(tmp1);
			p.eegIDs = tmp2(I);
		end
		eeg = ft_read_data(inName, 'chanindx', p.eegIDs); 
		eeg = reshape(eeg, [p.nEEGChans, p.rawNumPoints])';
		fprintf('EEG were read in!\n');
		clear Y I tmp0 tmp1 tmp2;
	else
		eeg = [];
		fprintf('no EEG read in!\n');
	end
	fs = [fs{:}, {'eeg'}];
	
	%% 5. HMV read
	if ~isempty(hmvID) && flagWetRun
		p.nHMVChans = length(hmvID);
		if isnumeric(hmvID)
			p.hmvIDs = hmvID;
			p.hmvLabels	= hdr.label(hmvID);
		elseif iscell(hmvID)
			p.hmvLabels	= hmvID;
			[tmp0, tmp1, tmp2] = intersect(hmvID, hdr.label);
			[Y, I] = sort(tmp1);
			p.hmvIDs = tmp2(I);
		end
		hmv = ft_read_data(inName, 'chanindx', p.hmvIDs); 
		hmv = reshape(hmv, [p.nHMVChans, p.rawNumPoints])';
		fprintf('head motion data were read in!\n');
		clear Y I tmp0 tmp1 tmp2;
		if strcmp(outName, 'hmv_only')	%this is a secret code!
			flagWetRun = false;
		end
	else
		hmv = [];
		fprintf('no head motion data read in!\n');
	end
	fs = [fs, {'hmv'}];
	
	%% 6. MEG read
	if flagWetRun
		meg = ft_read_data(inName, 'chanindx', p.megIDs); 
		sz = size(meg);
		if abs(p.numChannels - sz(1)) > 0.1 || ...
				abs(p.rawPointsPerTrial - sz(2)) > 0.1 || ...
				abs(p.rawNumTrials - sz(3)) > 0.1 || ...
				abs(p.rawNumPoints - prod(sz(2:end))) > 0.1
			error('Data size mismach with info!');
		end
		meg = reshape(meg, [p.numChannels, p.rawNumPoints])';
		fprintf('MEG were read in!\n');
	else
		meg = [];
		fprintf('no MEG read in!\n');
	end
	fs = [fs{:}, {'meg'}];
	
	%% 7. save data
	if flagWetRun > 0	
		if isempty(outName)
			[~, tmp, ~] = fileparts(inName);
			outName = [tmp, '_out.mat'];
			fprintf('outName not provided, assumed as %s!\n', outName);
		end
		try
			save(outName, fs{:});
		catch ME
			plt_myPrintLine('File Save Error!');
			fprintf('\n%s\n', ME.message);
		end
		fprintf('Data saved!\n');
	end
% 	ttt = clock;
% 	fprintf('Mission accomplished for %s @%02d:%02d:%02.0f...\n', ...
% 		inName, ttt(4:6));
	
end % end of function