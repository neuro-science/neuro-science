function [p, x, y] = mat_myFit (x0, y0)
% % % written 22/02/2018 by wp

	% % % start from polyfit 1
	if nargin < 3 %polyfit
		p = polyfit(x0, y0, 1);
		if nargout > 1
			m = min(x0);
			M = max(x0);
			mm = m - (M - m) * 0.1;
			MM = M + (M - m) * 0.1;
			x = linspace(mm, MM, 100);
			y = x * p(1) + p(2);
		end
	end
	

end
