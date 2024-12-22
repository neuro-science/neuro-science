function [t, p] = cmp12_ms2t(m1, m2, s1, s2, n1, n2)
% % % 17/07/14	written by wp: 
% % % mean and stand deviation of two samples
% % % number of cases

	% % % compute s
	if nargin < 5 || (n1 == n2)
		n = n1;
		s = sqrt((s1.^2 + s2.^2) / n);
		df = n - 1;
	else
		s = sqrt((((n1 - 1)*s1.^2 + (n2 - 1)*s2.^2)/(n1 + n2 - 2))...
			* (1/n1 + 1/n2));
		df = min(n1-1, n2-1);
	end
	
	% % % compute t	
	t = (m1 - m2) ./ s;
	
	% % % compute p if needed
	if nargout > 1
		p = 2 * tcdf(-abs(t), df);
	end
end % end of function