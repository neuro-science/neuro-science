function [artCorr, artSgnf] = meg_raw_corr_ICA_vs_EEG(ica, eeg)
% % %		25/08/14	re-written by wp
% % %		This function compute correlation of ica components to EOG and EKG
% % %		input:	ica(points, ics)
% % %					eeg(points, chs)


	%% 1. input check
	[npts, nICs] = size(ica);	% size of ica
	[npts2, nchs] = size(eeg);	% size of eeg
	if npts ~= npts2           % check consistence
		error('data size mismath between ICA and EEG!');
	end
	
	%% 2. compute
	% % % 	allocate space
	artCorr = zeros(nICs, nchs) + nan;
	artSgnf = zeros(nICs, nchs) + nan;
	
	% % % do it in a loop
	for ic = 1 : nICs
		for ie = 1 : nchs
			[rr, pp] = corrcoef(ica(:, ic), eeg(:, ie));
			artCorr(ic, ie) = rr(1, 2);	%[ica components, artfacts]
			artSgnf(ic, ie) = pp(1, 2);
		end
	end
end % end of functions