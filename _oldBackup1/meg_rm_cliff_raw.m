% % % This function is to remove cliff-like artifacts
% % % based on second order (differential) zscore.
% % % Input shall be: 
% % % data						: points x chans 
% % % infoTag (optional)	: char string of figure file prefix

function [data, out] = cmp13_rawmeg_rm_cliff(data, zThreshBase, zThreshRatio, zRange, nTorlerantRanges, infoTag)

% % % modified 23/07/2014	for more general use
% % % modified 4/9/2013	for more general use
	
	%% set parameters and prepare
	tic;
	% % % whether plot
	if nargin > 5 && ~isempty(infoTag)
		isPlot = 1;
		fid = fopen([infoTag, '.txt'], 'a');
	else
		isPlot = 0;
		fid = 1;
	end
	
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
	
	% % % data transform
	data = bsxfun(@minus, data, mean(data));	%remove channel mean
	
	% % show origin data plot
	if isPlot
		fh = figure('visible','off');
		plot(data);
		print(fh, '-dpng', [infoTag, '_allChans_before_Cliff_removal']);
		close(fh);
	end

	%% get bad points

	% % % get the change point
	z2 = zeros([npts, nchs]);
	z2(2:end, :) = abs(zscore((data(2:end, :) - data(1:end-1, :))));	%zscore of secondary change (derivative)
	theThresh = max(zThreshRatio * max(z2(:)), zThreshBase);	%conservative threshold
	[out.zPointsBad, out.zChansBad] = ind2sub([npts, nchs], find(z2 > theThresh));	%get indicis of bad points and channels
	clear z2;
	out.theChansBad = unique(out.zChansBad);	%get all bad channel ids
	out.nChansBad = length(out.theChansBad);	%number of bad channels
	out.nInterpolateDataPoints = 0; % the number of points that are interpolated

	% % % plot if needed
	if isPlot
		fh = figure('visible','off');
		plot(data(:, out.theChansBad));
		print(fh, '-dpng', [infoTag, '_selChans_before_Cliff_removal']);
		close(fh);
	end
	
	%% correction and interpolation
	out.nCliffs = out.nChansBad;

	% % % do the work channel by channel
	for iC = 1 : out.nChansBad
		myZChanBad = out.theChansBad(iC);	%the current bad channel
		myZPoints = out.zPointsBad(out.zChansBad == myZChanBad);	%get all above-thresh points for this channel
		nZPoints = length(myZPoints);	%number of above-thresh points

		% % % extend the range to cover sub-threshold artifacts 
		extPointsBad = repmat(myZPoints, [1 zRange*2 + 1]) + repmat(-zRange : zRange, [nZPoints, 1]);
		extPointsBad = unique(extPointsBad);	%remove overlappoint points
		extPointsBad(extPointsBad > npts) = npts;
		extPointsBad(extPointsBad < 1) = 1;
		
		dPB = ones(size(extPointsBad));	%we need derivatives
		dPB(2:end) = extPointsBad(2:end) - extPointsBad(1:end-1);	% get the derivatives of the points
		idxTurnPoints = find(dPB > nTorlerantRanges * zRange);	%get turning points index after large gaps
		theSegments(:, 1) = extPointsBad([1; idxTurnPoints]);	%start points of bad data segments
		theSegments(:, 2) = extPointsBad([idxTurnPoints - 1; end]);	%end points of bad data segments
		nSegments = size(theSegments, 1);	%number of segments
		% % % 	correction for data before first bad points
		thisPiece = data(1 : theSegments(1, 1), myZChanBad);	%take the piece of data
		data(1 : theSegments(1, 1), myZChanBad) = thisPiece - mean(thisPiece);	%rereference of new mean
		clear thisPiece;
		% % % 	correction for data after last bad points
		thisPiece = data(theSegments(end, 2) : end, myZChanBad);	%take the piece of data
		data(theSegments(end, 2) : end, myZChanBad) = thisPiece - mean(thisPiece);	%rereference of new mean
		clear thisPiece;
		% % % 	correction for data within gaps of bad data
		if nSegments > 1
			for iS = 2 : nSegments
				thisPiece = data((theSegments(iS - 1, 2) : theSegments(iS, 1)), myZChanBad);
				data((theSegments(iS - 1, 2) : theSegments(iS, 1)), myZChanBad) = thisPiece - mean(thisPiece);
				clear thisPiece;
			end
			out.nCliffs = out.nCliffs + nSegments - 1;
		end
		% % % 	interpolation for data within the bad range
		for iS = 1 : nSegments
			thisPiece = data((theSegments(iS, 1) : theSegments(iS, 2)), myZChanBad);
			data((theSegments(iS, 1) : theSegments(iS, 2)), myZChanBad) = ...
				linspace(thisPiece(1), thisPiece(end), length(thisPiece));	%linear interpolation
			out.nInterpolateDataPoints = out.nInterpolateDataPoints + length(thisPiece);
			clear thisPiece;
		end
		% % % save for further use
		out.Segment{iC} = theSegments;
		clear theSegments myZChanBad myZPoints nZPoints extPointsBad dPB idxTurnPoints nSegments;
	end

	%% return and clean
	% % % output message
	xpts = unique(out.zPointsBad);
	out.rCliffs = 100 * length(xpts) ./ npts;
	fprintf(fid, 'After %4.2f minutes, %d cliffs were found and corrected, %4.3f%% data were interpolated,\n', ...
		toc/60, out.nCliffs, out.rCliffs);
	clear idx tmp xpts;

	% % % show corrected data plot
	if isPlot
		fh = figure('visible','off');
		plot(data);
		print(fh, '-dpng', [infoTag, '_allChans_after_Cliff_removal']);
		clf;
		plot(data(:, out.theChansBad));
		print(fh, '-dpng', [infoTag, '_selChans_after_Cliff_removal']);
		close(fh);
	end
	clear minIndBad maxIndBad indRangeFull preBadData posBadData;
	if isPlot
		fclose(fid);
	end
end % end of function