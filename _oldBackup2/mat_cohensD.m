
function D = mat_cohensD(d1, d2, pairFlag)

% % % updated on 18/07/2019 by wp, for multi-dimensional data
% % % written on 13/02/2019 by wp, compute cohen's D for two dataset
% % % fist dimension for cases

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
	elseif numel(size(d1)) == 2 && (size(d1, 2) == 1 || size(d1, 1) == 1) && numel(size(d2)) == 2 && (size(d2, 2) == 1 || size(d2, 1) == 1)
		m1 = mean(d1);
		m2 = mean(d2);
		s1 = std(d1);
		s2 = std(d2);
		n1 = numel(d1);
		n2 = numel(d2);
		s = sqrt(((n1-1)*s1.^2 + (n2-1)*s2.^2)/(n1 + n2 - 2));
		D = (m1 - m2) / s;
	elseif all(size(d1) == size(d2))
		m1 = mean(d1);
		m2 = mean(d2);
		s1 = std(d1);
		s2 = std(d2);
		n1 = size(d1, 1);
		s = sqrt(((n1-1)*s1.^2 + (n1-1)*s2.^2)/(n1*2 - 2));
		D = (m1 - m2) ./ s;
	else
		error('data size mismatch!');
	end

end %end of function