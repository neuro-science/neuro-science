% % % This function is to remove between block noise
% % % Input shall be: 
% % % data							: points x chans 
% % % theSegs						: nSegments x 2 (start, end)
% % % tolStart						: tolerance data points for start
% % % tolEnd						: tolerance data points for end
% % % updated 03/10/2023 by wp
% % % written 11/02/2022 by wp

function [data, out] = meg_raw_rm_breaks (data, theSegs, tolStart, tolEnd, flagMove)
	
	%% 1. input check and headers
	% % % too many or too few inputs?
	if nargin < 2
		error('Data (points x chans) and theSegs (n x 2) are needed!');
	elseif nargin < 3 || isempty(tolStart)
		tolStart = 2400;
		tolEnd = 2400;
	elseif nargin < 4 || isempty(tolEnd)
		tolEnd = tolStart;
	elseif nargin < 5 || isempty(flagMove)
		flagMove = true;
	elseif nargin > 5
		error('Too many inputs!');
	end
	
	% % % data size check
	[npts, nchs] = size(data);
	if npts < nchs
		fprintf('\nWarning: I assume input data is "chans x points".\n');
		fprintf('I will change it to "points x chans"...\n');
		data = data';
		[npts, nchs] = size(data);
		fprintf('Done.\n\n');
	end		
	
	% % % check the theSegs size
	sz = size(theSegs);
	if abs(sz(2) - 2) > 0.1
		if abs(sz(1) - 2) < 0.1
			fprintf('\nWarning: I assume input data is [start; end] x n.\n');
			fprintf('I will change it to the default: n x [start, end]...\n');
			tmp = theSegs;
			theSegs = theSegs';
			sz = size(theSegs);
		else
			error('The size of theSegs should be n x [start, end]!');
		end
	else
		tmp = theSegs';
	end
	
	% % % check the validity of theSegs	
	tmp1 = tmp(:);
	tmp2 = sort(tmp1, 'ascend');
% 	tmp3 = any(theSegs(:, 2) < theSegs(:, 1));
	if ~isequal(tmp1, tmp2)	|| tmp1(1) < 1 || tmp1(end) > npts %|| tmp3
		error('the segements were implausible!')
	else
		clear tmp tmp1 tmp2;	% tmp3
	end
	
	if flagMove
		%% 2. do the removal in move mode
		% % % remove the break for the first segment	
		if (theSegs(1, 1) - tolStart - 1) < 0.1 % First point is 1
			data = bsxfun(@minus, data, data(1, :));
			out = [];
		else
			data(theSegs(1, 1) - tolStart : end, :) = bsxfun(@minus, ...
				data(theSegs(1, 1) - tolStart : end, :), data(theSegs(1, 1) - tolStart, :));

			% % % piece before first block		
			out = [1, theSegs(1, 1) - 1];
			data(1 : theSegs(1, 1) - tolStart - 1, :) = 0;
		end

		% % % work on the break for the other segments	
		for ib = 2 : sz(1)
			% % % fill breaks	
			data(theSegs(ib - 1, 2) + tolEnd + 1 : theSegs(ib, 1) - tolStart - 1, :) = ...
				repmat(data(theSegs(ib - 1, 2) + tolEnd, :), ...
				[theSegs(ib, 1) - theSegs(ib - 1, 2) - tolStart - tolEnd - 1, 1]);

			% % % difference for start point and previous end
			tmpDiff = data(theSegs(ib, 1) - tolStart, :) - ...
				data(theSegs(ib - 1, 2) + tolEnd, :);

			% % % remove above difference for this segment
			data(theSegs(ib, 1) - tolStart : end, :) = bsxfun(@minus, ...
				data(theSegs(ib, 1) - tolStart : end, :), tmpDiff);

			% % % offer the inter-block points
			out = [out; theSegs(ib - 1, 2) + 1, theSegs(ib, 1) - 1];
		end

		% % % part after the last block
		if (theSegs(sz(1), 2) - npts) < -0.1
			out = [out; theSegs(sz(1), 2) + 1, npts];
		end
		if theSegs(sz(1), 2) + tolEnd < npts
			data(theSegs(sz(1), 2) + tolEnd + 1 : npts, :) = repmat(...
				data(theSegs(sz(1), 2) + tolEnd, :), [npts - theSegs(sz(1), 2) - tolEnd, 1]);
		end
	else
		%% 3. do the removal in retain mode
		% % % remove the break for the first segment	
		if (theSegs(1, 1) - tolStart - 1) < 0.1 % First point is 1
			out = [];
		else
			data(1 : theSegs(1, 1) - tolStart - 1, :) = repmat( ...
				data(theSegs(1, 1) - tolStart, :), [theSegs(1, 1) - tolStart - 1, 1]);

			% % % piece before first block		
			out = [1, theSegs(1, 1) - 1];
		end

		% % % work on the break for the other segments	
		for ib = 2 : sz(1)
			% % % fill breaks	with linear interploration
			N = theSegs(ib, 1) - theSegs(ib - 1, 2) - tolStart - tolEnd;
			S = repmat((0 : N) / N, [nchs, 1])';
			D = bsxfun(@mtimes, data(theSegs(ib, 1) - tolStart, :) - ...
				data(theSegs(ib - 1, 2) + tolEnd, :), S);
			data(theSegs(ib - 1, 2) + tolEnd : theSegs(ib, 1) - tolStart, :) = ...
				bsxfun(@plus, data(theSegs(ib - 1, 2) + tolEnd, :), D);
			clear N S D;
			% % % offer the inter-block points
			out = [out; theSegs(ib - 1, 2) + 1, theSegs(ib, 1) - 1];
		end

		% % % part after the last block
		if (theSegs(sz(1), 2) - npts) < -0.1
			out = [out; theSegs(sz(1), 2) + 1, npts];
		end
		if theSegs(sz(1), 2) + tolEnd < npts
			data(theSegs(sz(1), 2) + tolEnd + 1 : npts, :) = repmat(...
				data(theSegs(sz(1), 2) + tolEnd, :), [npts - theSegs(sz(1), 2) - tolEnd, 1]);
		end
	end
end % end of functiontheSegs(ib, 1) - tolStart - theSegs(ib - 1, 2) - tolEnd - 1