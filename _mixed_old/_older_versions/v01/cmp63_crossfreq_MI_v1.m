function [pac, pac_r] = cmp63_crossfreq_MI(a, A, N, M, MM)
% % % 12/01/17	written by wp 
% % % compute modulation Index based on entropy measures
% % % 	a: angle of low frequency data (points x channels)
% % % 	A: Amplitude of high frequency data (points x channels)
% % %		N: number of bins within a cycle, default 12

	%% 1. check inputs
	% % % number of bins
	if nargin > 3
		randFlag = true;
		if nargin < 5
			MM = 100;
		end
	else
		randFlag = false;
		pac_r = [];
	end
	if nargin < 3
		N = 12;
	end
	% % % input data
	if nargin < 2
		error('We need at least two inputs: angle and Amplitude!');
	end
	% % % data size
	s1 = size(a);
	s2 = size(A);
	if abs(s1(1) - s2(1)) > 0.1 || length(s1) > 2 || length(s2) > 2
		error('data size mismatch!');
	end
	pac = zeros(s1(2), s2(2)) + nan;
	
	%% 2. work on the data
	borders = linspace(-pi, pi, N+1);
	if randFlag % with permutation
		% % % initialize
		pac_r = zeros([s1(2), s2(2), MM]) + nan;
		ids = zeros(s1(1), MM, 'uint32');
		nTrials = s1(1) / M;	
		for m = 1 : MM
			tmp1 = randperm(nTrials);
			tmp2 = bsxfun(@plus, tmp1 * M - M, (1 : M)');
			ids(:, m) = tmp2(:);
		end
		clear tmp1 tmp2;
		for c1 = 1 : s1(2)
			[x1, x2] = histc(a(:, c1), borders);
			d1 = zeros(N, s2(2));
			parfor ib = 1 : N
				d1(ib, :) = mean(A(x2 == ib, :), 1);
			end
			pj = bsxfun(@rdivide, d1, sum(d1, 1));
			pac(c1, :) = 1 + sum(pj .* log(pj)) ./ log(N);
			% % % do statistics for shuffled data
			for m = 1 : MM
				A1 = A(ids(:, m), :);
				d1 = zeros(N, s2(2));
				parfor ib = 1 : N
					d1(ib, :) = mean(A1(x2 == ib, :), 1);
				end
				pj = bsxfun(@rdivide, d1, sum(d1, 1));
				pac_r(c1, :, m) = 1 + sum(pj .* log(pj)) ./ log(N);
			end
		end
	else % without permutation
		for c1 = 1 : s1(2)
			[x1, x2] = histc(a(:, c1), borders);
			d1 = zeros(N, s2(2));
			parfor ib = 1 : N
				d1(ib, :) = mean(A(x2 == ib, :), 1);
			end
			pj = bsxfun(@rdivide, d1, sum(d1, 1));
			pac(c1, :) = 1 + sum(pj .* log(pj)) ./ log(N);
		end
	end		
end %end of function

