function cmp61_raw_proTimingTestMEG (inName, outName, eegID, rangeT)
	% % % This function will deal with the trigger timing test 
	% % % written by peng wang @09/08/2016

	%% 1. read the data and pre-processing
	% % % load data and save
	[p, hdr, evt] = fcv02_ds2mat(inName, outName, eegID);
	load(outName);
	save(outName, 'eeg', 'evt', 'hdr');
	
	% % % filter the data
	[tmp1, tmp2] = cmp10_rawmeg_filter(eeg, hdr.Fs, 55, 45, 2);
	eeg2 = zscore(tmp1);
	x = (1 : length(eeg2))./hdr.Fs; %time in seconds
	[npts, nchs] = size(eeg2);
	% % % get the triggers
	ct = 0;
	for k = 1 : length(evt)
		if strcmp(evt(k).type, 'UPPT001')
			ct = ct + 1;
			trg(ct, 1) = evt(k).sample;
			trg(ct, 2) = evt(k).value;
		end
	end
	
	% % % epoch the data 
	rangeS = round(rangeT * hdr.Fs);
	xr = linspace(rangeT(1), rangeT(2), rangeS(2) - rangeS(1) + 1);
	ntrs = size(trg, 1);
	eeg3 = zeros(rangeS(2) - rangeS(1) + 1, nchs, ntrs) + nan;
	for it = 1 : ntrs
		if trg(it, 1) + rangeS(1) < 1
			eeg3(2 - trg(it, 1) - rangeS(1) : end, :, it) = eeg2(1 : trg(it, 1) + rangeS(2), :);
		elseif trg(it, 1) + rangeS(2) > npts
			eeg3(1 : npts - trg(it, 1) - rangeS(1) + 1, :, it) = eeg2(trg(it, 1) + rangeS(1) : end, :);
		else
			eeg3(:, :, it) = eeg2(trg(it, 1) + rangeS(1) : trg(it, 1) + rangeS(2), :);
		end
	end
	
% 	% % % shift the data
% 	for k = 1 : 6
% 		ea(:, :, k) = [zeros(k, 2); eeg2(1 : end - k, :)];
% 		eb(:, :, k) = [eeg2(k + 1 : end, :); zeros(k, 2)];
% 	end
% 	ec = zscore((eb(:, :, 6) - ea(:, :, 6))/12);
% 	
% 	% % % find change point?	
% 	ea = [zeros(1, 2); ec(1 : end - 1, :)];
% 	eb = [ec(2 : end, :); zeros(1, 2)];
% 	ee = ec > ea & ec > eb & ec > th;
% 	es = ec < ea & ec < eb & ec < -th;
% 	ee = find(ee(:, 2)) ./ hdr.Fs; 
% 	es = find(es(:, 2)) ./ hdr.Fs;
	
	% % % clean up
	save(outName, 'eeg2', 'eeg3', 'ntrs', 'nchs', 'trg', 'x', 'xr', '-append');
% 	save(outName, 'eeg2', 'ee', 'es', 'ec', 'trg', 'x', '-append');
	clear eeg evt hdr meg p inName outName ea eb eeg2 eeg3;
end


