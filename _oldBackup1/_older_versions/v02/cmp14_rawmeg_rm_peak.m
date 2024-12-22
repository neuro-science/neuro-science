% % % This function is to remove peak like artifact based on z-score
% % % Input and output data are both continuous
function [data, out] = cmp14_rawmeg_rm_peak(data, myThresh, myZRange, infoTag)

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
		fprintf(fid, fid, '\nWarning: I assume your input data is "chans x points".\n');
		fprintf(fid, fid, 'I will change it to "points x chans"...\n');
		data = data';
		[npts, nchs] = size(data);
		fprintf(fid, fid, 'Done.\n\n');
	end		

	%% find peak and label as not usable
	z = abs(zscore(data));	%zscore was done in all channels, data shall be continuous
	[out.pts, out.chs] = ind2sub([npts, nchs], find(z > out.zThresh));	%get the indices

	ptsBad = [];
	if ~isempty (out.chs)	% peaks detected

		% % % interpolate for filtering, the data will not be used anyway
		ptsBad = unique(out.pts);	%remove overlappoint points
		dPB = ones(size(ptsBad));	%we need derivatives
		dPB(2:end) = ptsBad(2:end) - ptsBad(1:end-1);	% get the derivatives of the points
		idxTurnPoints = find(dPB > out.zRange);	%get turning points index after large gaps
		out.theSegments(:, 1) = ptsBad([1; idxTurnPoints]);	%start points of bad data segments
		out.theSegments(:, 2) = ptsBad([idxTurnPoints - 1; end]);	%end points of bad data segments
		nSegments = size(out.theSegments, 1);	%number of segments
		out.theChansBad = unique(out.chs);

		% % % 	correction for data before first bad points - remove mean only
		thisPiece = data(1 : out.theSegments(1, 1), out.theChansBad);	%take the piece of data
		data(1 : out.theSegments(1, 1), out.theChansBad) = ...
			bsxfun(@minus, thisPiece, mean(thisPiece));	%rereference of new mean
		clear thisPiece;

		% % % 	correction for data after last bad points - remove mean only
		thisPiece = data(out.theSegments(end, 2) : end, out.theChansBad);	%take the piece of data
		data(out.theSegments(end, 2) : end, out.theChansBad) = ...
			bsxfun(@minus, thisPiece, mean(thisPiece));	%rereference of new mean
		clear thisPiece;

		% % % 	correction for data within gaps between bad data - remove mean only
		if nSegments > 1
			for iS = 2 : nSegments
				thisPiece = data((out.theSegments(iS - 1, 2) : out.theSegments(iS, 1)), out.theChansBad);
				data((out.theSegments(iS - 1, 2) : out.theSegments(iS, 1)), out.theChansBad) = ...
					bsxfun(@minus, thisPiece, mean(thisPiece));	%rereference of new mean
				clear thisPiece;
			end
		end
		
		% % % 	interpolation for data within the bad range
		for iS = 1 : nSegments
			N = out.theSegments(iS, 2) - out.theSegments(iS, 1) + 1;	%length of data
			if out.theSegments(iS, 1) > 1
				m = data(out.theSegments(iS, 1) - 1, out.theChansBad);		%start data, one point before
			else
				m = mean(data(:, out.theChansBad));	%take mean if this is the first point
			end
			if out.theSegments(iS, 2) < size(data, 1)
				M = data(out.theSegments(iS, 2) + 1, out.theChansBad);		%end data, one point after
			else
				M = mean(data(:, out.theChansBad));	%take mean if it is the last point
			end
			% % % 	linear interpolation
			data((out.theSegments(iS, 1) : out.theSegments(iS, 2)), out.theChansBad) = ...
				bsxfun(@plus, m, ((1 : N) / (N + 1))' * (M - m));	%not the direction of vectors before *
			clear M m N;
		end
	else
		fprintf(fid, '===No peaks detected in this dataset!===\n');
	end
	
	%% return and clean
	% % % show corrected data plot
	if isPlot
		fh = figure('visible','off');
		plot(data);
		print(fh, '-dpng', [infoTag, '_allChans_after_Peak_removal']);
		close(fh);
	end
	out.rjRate = 100 * length(ptsBad) ./ npts;
	
	fprintf(fid, 'After %4.2f minutes, %5.2f%% data was labeled as rejected because of jump artifacts!\n', toc/60, out.rjRate);
	if isPlot
		fclose(fid);
	end
end