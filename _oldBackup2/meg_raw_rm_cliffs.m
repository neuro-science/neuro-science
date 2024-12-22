% % % This function is to remove cliff-like artifacts
% % % based on second order (differential) zscore.
% % % Input shall be: 
% % % data						: points x chans 
% % % infoTag (optional)	: char string of figure file prefix

% % % rewriten 04/10/2023	for break algorithm
% % % modified 23/07/2014	for more general use
% % % modified 04/09/2013	for more general use

function [data, out] = meg_raw_rm_cliff(data, ...
	zThreshBase, zThreshRatio, zRange, nTorlerantRanges)
	
	%% 1. set parameters and prepare
	
	% % % data size
	[npts, nchs] = size(data);
	if npts < nchs
		fprintf(fid, '\nWarning: I assume your input data is "chans x points".\n');
		fprintf(fid, 'I will change it to "points x chans"...\n');
		data = data';
		[npts, nchs] = size(data);
		fprintf(fid, 'Done.\n\n');
	end		

	% % % other parameters	
	if nargin < 5 || isempty(nTorlerantRanges)
		nTorlerantRanges = 5;
	end
	if nargin < 4 || isempty(zRange)
		zRange = 30;
	end
	if nargin < 3 || isempty(zThreshRatio)
		zThreshRatio = 0.1;
	end
	if nargin < 2 || isempty(zThreshBase)
		zThreshBase = 100;
	end

	%% 2. get bad points

	% % % get the change point
	z2 = zeros([npts, nchs]);
	% % % zscore of secondary change (derivative)
	z2(2:end, :) = abs(zscore((data(2:end, :) - data(1:end-1, :))));	
	% % % conservative threshold
	out.theThresh = max(zThreshRatio * max(z2(:)), zThreshBase);	
	% % % get indicis of bad points and channels
	out.zPointsBad = find(any(z2 > out.theThresh, 2));
	% % % quit if there was no cliff	
	if isempty(out.zPointsBad)
		fprintf('No cliffs were found!\n');
		return;
	end
	clear z2;
	
	%% 3. correction and interpolation

	% % % extended bad points to cover ranges
	extPointsBad = bsxfun(@plus, out.zPointsBad, -zRange : zRange);
	extPointsBad(extPointsBad > npts) = npts;
	extPointsBad(extPointsBad < 1) = 1;
	extPointsBad = unique(extPointsBad);	

	% % % break points into pieces	
	% % % we need derivatives
	dPB = ones(size(extPointsBad));	
	% % % get the derivatives of the points
	dPB(2:end) = extPointsBad(2:end) - extPointsBad(1:end-1);	
	% % % get turning points index after large gaps
	idxTurnPoints = find(dPB > nTorlerantRanges * zRange);	
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
	
	% % % remove the broken pieces
	[data, out.theSegments] = meg_raw_rm_breaks (data, theSegments, 0, 0);
	
	%% 4. outputs
	% % % number of Cliffs
	out.nCliffs = size(out.theSegments, 1);
	% % % number of bad data points
	out.nInterpolateDataPoints = sum(out.theSegments(:, 2) - ...
		out.theSegments(:, 1)) + out.nCliffs;
	% % % ratio of bad data
	out.rCliffs = 100 * out.nInterpolateDataPoints ./ npts;
	% % % output message
	fprintf('%d cliffs (%4.3f%%) were found and corrected!\n', ...
		out.nCliffs, out.rCliffs);
	
end % end of function