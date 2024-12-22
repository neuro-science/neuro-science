	function t = mat_tValue4Arrays(data, dim, sqFlag)
	% % % 22/06/17	modified by wp: rename for general use
	% % % 26/05/14	modified by wp: use nan for the missing values
	% % % 12/05/14 written by wp

		% % % check paras
		if nargin < 3 || isempty(sqFlag)
			sqFlag = 1;
		end
		if nargin < 2 || isempty(dim) %default last dim
			dim = length(size(data));
		end

		% % % compute
		n = size(data, dim);
		m = nanmean(data, dim);
		s = nanstd(data, 0, dim);
		t = sqrt(n) * (m./s);

		% % % squeeze if needed
		if sqFlag
			t = squeeze(t);
		end
	end
