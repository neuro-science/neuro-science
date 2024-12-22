% % % This function is to remove peak like artifact based on z-score
% % % Input and output data are both continuous

% % % rewritten 04/10/2023 by wp for break algorithm
% % % rewritten 23/07/2014 by wp for more general purpose
function [data, out] = meg_raw_rm_peaks(data, myThresh, myZRange, ...
	flagDry,	nTorlerantRanges)

	%% 1. initialization

	% % % check parameters
	if nargin < 5 || isempty(nTorlerantRanges)
		nTorlerantRanges = 5;	% The threshold was derived from previous manual work
	end
	
	if nargin < 4 || isempty(flagDry)
		flagDry = false;	% The threshold was derived from previous manual work
	end
	
	if nargin < 3 || isempty(myZRange)
		myZRange = 600;	% The threshold was derived from previous manual work
	end
	
	if nargin < 2 || isempty(myThresh)
		myThresh = 6;	% The threshold was derived from previous manual work
	end
	
	% % % data size
	[npts, nchs] = size(data);
	if npts < nchs
		fprintf(fid, fid, '\nWarning: I assume your input data is "chans x points".\n');
		fprintf(fid, fid, 'I will change it to "points x chans"...\n');
		data = data';
		[npts, nchs] = size(data);
		fprintf(fid, fid, 'Done.\n\n');
	end		

	%% 2. find peak and label as not usable
	% % % zscore was done in all channels, data shall be continuous
	z = abs(zscore(data));	
	% % % get indicis of bad points and channels
	out.zPointsBad = find(any(z > myThresh, 2));
	% % % quit if there was no cliff	
	if isempty(out.zPointsBad)
		fprintf('No peaks were found!\n');
		return;
	end
	
	%% 3. correction and interpolation

	% % % extended bad points to cover ranges
	extPointsBad = bsxfun(@plus, out.zPointsBad, -myZRange : myZRange);
	extPointsBad(extPointsBad > npts) = npts;
	extPointsBad(extPointsBad < 1) = 1;
	extPointsBad = unique(extPointsBad);	

	% % % break points into pieces	
	% % % we need derivatives
	dPB = ones(size(extPointsBad));	
	% % % get the derivatives of the points
	dPB(2:end) = extPointsBad(2:end) - extPointsBad(1:end-1);	
	% % % get turning points index after large gaps
	idxTurnPoints = find(dPB > nTorlerantRanges * myZRange);	
	% % % start points of good data segments
	s1 = extPointsBad(idxTurnPoints - 1) + 1;
	% % % end points of good data segments
	s2 = extPointsBad(idxTurnPoints) - 1;
	% % % start NOT with point 1
	if abs(extPointsBad(1) - 1) > 0.1	
		s1 = [1; s1];
		s2 = [extPointsBad(1) - 1; s2];
	end
	% % % end NOT with point end
	if abs(extPointsBad(end) - npts) > 0.1	
		s1 = [s1; extPointsBad(end) + 1];
		s2 = [s2; npts];
	end
	% % % end points of bad data segments
	theSegments = [s1, s2];	
	% % % clear	up
	clear dPB extPointsBad idxTurnPoints s1 s2;

	% % % return number of segments for dry run
	if flagDry
		out.nPeaks = size(theSegments, 1);
		fprintf('%d peaks offered as a dry run!\n', out.nPeaks);
		return;
	end
	
	% % % remove the broken pieces in retain mode
	[data, out.theSegments] = meg_raw_rm_breaks (data, theSegments, 0, 0, false);
	
	
	%% 4. outputs
	% % % number of Cliffs
	out.nPeaks = size(out.theSegments, 1);
	% % % number of bad data points
	out.nInterpolateDataPoints = sum(out.theSegments(:, 2) - out.theSegments(:, 1)) ...
		+ out.nPeaks;
	% % % ratio of bad data
	out.rjRate = 100 * out.nInterpolateDataPoints ./ npts;
	% % % output message
	fprintf('%d peaks (%4.3f%%) were found and corrected!\n', ...
		out.nPeaks, out.rjRate);
end