% % % This function is to remove peak like artifact based on z-score
% % % Input and output data are both continuous
function [data, out] = meg_raw_rm_peak(data, myThresh, myZRange, infoTag)

% % % rewritten 16/09/2022 by wp compensation with average
% % % rewritten 12/03/2022 by wp compensation as breaks (neighbor)
% % % rewritten 23/07/2014 by wp for more general purpose

	%% initialization
	tic;
	% % % check whether log needed	
	if nargin > 3 && ~isempty(infoTag)
		isPlot = 1;
		fid = fopen([infoTag, '.txt'], 'a');
	else
		isPlot = 0;
		fid = 1;
	end

	% % % check parameters
	if nargin < 3 || isempty(myZRange)
		out.zRange = 600;	% The threshold was derived from previous manual work
	else
		out.zRange = myZRange;
	end
	
	if nargin < 2 || isempty(myThresh)
		out.zThresh = 6;	% The threshold was derived from previous manual work
	else
		out.zThresh = myThresh;
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

	%% find peak and label as not usable
	z = abs(zscore(data));	%zscore was done in all channels, data shall be continuous
	[out.pts, out.chs] = ind2sub([npts, nchs], find(z > out.zThresh));	%get the indices
	out.theChansBad = unique(out.chs);

	if ~isempty (out.theChansBad)	% peaks detected

		% % % interpolate for filtering, the data will not be used anyway
		ptsBad = unique(out.pts);	%remove overlappoint points
		extPointsBad = bsxfun(@plus, ptsBad, -out.zRange : out.zRange);
		extPointsBad(extPointsBad > npts) = npts;
		extPointsBad(extPointsBad < 1) = 1;
		extPointsBad = unique(extPointsBad);	%remove overlappoint points
		out.rjRate = 100 * numel(extPointsBad) ./ npts;
		
		% % % do cleanning average outside of bad points		
		data2 = data;
		data2(extPointsBad, :) = [];
		m = nanmean(data2, 1);
		clear data2;
		data(extPointsBad, :) = repmat(m, [numel(extPointsBad) 1]);
		
		dPB = ones(size(extPointsBad));	%we need derivatives
		dPB(2:end) = extPointsBad(2:end) - extPointsBad(1:end-1);	% get the derivatives of the points
		idxTurnPoints = find(dPB > 1);	%get turning points index after 1 gaps
		out.theSegments(:, 1) = extPointsBad([1; idxTurnPoints]);	%start points of bad data segments
		out.theSegments(:, 2) = extPointsBad([idxTurnPoints - 1; end]);	%end points of bad data segments
		clear m ptsBad extPointsBad idxTurnPoints dPB;
		
	else
		out.rjRate = 0;
	end
	%% return and clean
	% % % show corrected data plot
	if isPlot
		fh = figure('visible','off');
		plot(data);
		print(fh, '-dpng', [infoTag, '_allChans_after_Peak_removal']);
		close(fh);
		fclose(fid);
	end
	
	fprintf(fid, 'After %4.2f minutes, %5.2f%% data was labeled as rejected because of jump artifacts!\n', toc/60, out.rjRate);
end