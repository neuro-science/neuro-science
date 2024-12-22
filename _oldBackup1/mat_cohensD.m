function D = mat_cohensD(d1, d2, pairFlag)
% % % written on 13/02/2019 by wp, compute cohen's D for two dataset

	%% initialize
	if nargin < 3 || isempty(pairFlag)
		pairFlag = false;
	end
	
	if nargin < 2
		pairFlag = true;
	elseif pairFlag
		d1 = d1 - d2;
	end
	
	if pairFlag
		m = mean(d1);
		s = std(d1);
		D = m / s;
	else
		m1 = mean(d1);
		m2 = mean(d2);
		s1 = std(d1);
		s2 = std(d2);
		n1 = numel(d1);
		n2 = numel(d2);
		s = sqrt(((n1-1)*s1.^2 + (n2-1)*s2.^2)/(n1 + n2 - 2));
		D = (m1 - m2) / s;
	end

end %end of function