% % % This function is to remove cliff-like artifacts
% % % based on second order (differential) zscore.
% % % Input shall be: 
% % % data						: points x chans 
% % % infoTag (optional)	: char string of figure file prefix

function [data, out] = meg_raw_rm_cliff(data, zThreshBase, zThreshRatio, zRange, infoTag)

% % % modified 11/02/2022	for better performance
% % % modified 23/07/2014	for more general use
% % % modified 04/09/2013	for more general use
	
	%% 1. set parameters and prepare
	tic;
	% % % whether plot
	if nargin > 4 && ~isempty(infoTag)
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
	if nargin < 4 || isempty(zRange)
		zRange = 30;
	end
	if nargin < 3 || isempty(zThreshRatio)
		zThreshRatio = 0.1;
	end
	if nargin < 2 || isempty(zThreshBase)
		zThreshBase = 100;
	end
	
	% % show origin data plot
	if isPlot
		fh = figure('visible','off');
		plot(data);
		print(fh, '-dpng', [infoTag, '_allChans_before_Cliff_removal']);
		close(fh);
	end

	%% 2. get bad points
	% % % get the change point
	z2 = zeros([npts, nchs]);
	z2(2:end, :) = abs(zscore((data(2:end, :) - data(1:end-1, :))));	%zscore of secondary change (derivative)
	theThresh = max(zThreshRatio * max(z2(:)), zThreshBase);	%conservative threshold
	[out.zPointsBad, out.zChansBad] = ind2sub([npts, nchs], find(z2 > theThresh));	%get indicis of bad points and channels
	clear z2 theThresh;
	
	if isempty(out.zPointsBad)
		% % % output message
		fprintf(fid, 'After %4.2f minutes, NO cliffs were found.\n', toc/60);
	else
		out.theChansBad = unique(out.zChansBad);	%get all bad channel ids
		out.nChansBad = numel(out.theChansBad);	%number of bad channels
		% % % extend the range to cover sub-threshold artifacts 
		myZPoints = unique(out.zPointsBad);	%get all bad data points %[nZPoints 1]
		extPointsBad = bsxfun(@plus, myZPoints, -zRange : zRange);
		extPointsBad(extPointsBad > npts) = npts;
		extPointsBad(extPointsBad < 1) = 1;
		extPointsBad = unique(extPointsBad);	%remove overlappoint points
		out.rCliffs = 100 * numel(extPointsBad) ./ npts;

		% % % plot if needed
		if isPlot
			fh = figure('visible','off');
			plot(data(:, out.theChansBad));
			print(fh, '-dpng', [infoTag, '_selChans_before_Cliff_removal']);
			close(fh);
		end

		%% correction and interpolation with all channels % updated in newer version 2022
		% % % find the turning points	
		dPB = ones(size(extPointsBad));	%we need derivatives
		dPB(2:end) = extPointsBad(2:end) - extPointsBad(1:end-1);	% get the derivatives of the points
		idxTurnPoints = find(dPB > 1);	%get turning points index after 1 gap
		% % % find the segments
		out.theSegments(:, 1) = extPointsBad([1; idxTurnPoints]);	%start points of bad data segments
		out.theSegments(:, 2) = extPointsBad([idxTurnPoints - 1; end]);	%end points of bad data segments
		out.nSegments = size(out.theSegments, 1);	%number of segments
		clear dPB idxTurnPoints extPointsBad myZPoints;
		% % % do cleaning
		data = meg_raw_rm_breaks(data, out.theSegments, out.theChansBad);
		data = bsxfun(@minus, data, mean(data));

		% % % show corrected data plot
		if isPlot
			fh = figure('visible','off');
			plot(data);
			print(fh, '-dpng', [infoTag, '_allChans_after_Cliff_removal']);
			clf;
			plot(data(:, out.theChansBad));
			print(fh, '-dpng', [infoTag, '_selChans_after_Cliff_removal']);
			close(fh);
			fclose(fid);
		end

		% % % output message
		fprintf(fid, 'After %4.2f minutes, %d cliffs were found and corrected, %4.3f%% data were interpolated,\n', ...
			toc/60, out.nSegments, out.rCliffs);
	end
end % end of function