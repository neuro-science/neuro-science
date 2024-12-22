function [pac, pac_r] = cmp63_crossfreq_MI(a, A, N, M, MM)
% % % 16/01/17	modified by wp, increase speed at price of memory 
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
	
	%% 2. work on the data
	% % % work on the phase channel
	borders = linspace(-pi, pi, N+1);
	da = zeros(size(a)) + nan;
	for c1 = 1 : s1(2)
		[tmp1, da(:, c1)] = histc(a(:, c1), borders);
	end
	da = da';
	
	% % % logical sets
	X = false([size(a), N]);
	for ib = 1 : N
		X(:, :, ib) = da == ib;
	end
	
	% % % take the amp channel in
	dA = ~isnan(A);
	A(dA) = 0;
	
	% % % mean distribution
	sn = bsxfun(@mtimes, X, double(dA));
	pj = zeros(s1(2), s2(2), N) + nan;
	for ib = 1 : N
		X = da == ib;
		sn = X * double(dA);
		sd = X * A;
		pj(:, :, ib) = sd ./ sn;
	end
	pm = bsxfun(@rdivide, pj, sum(pj, 3));
	
	% % % summary
	pac = 1 + sum(pm .* log(pm), 3) ./ log(N);
	clear pj pm sn sd;
	
	%% 3. statistics 
	if randFlag % with permutation
		% % % initialize
		pac_r = zeros([s1(2), s2(2), MM]) + nan;
		nTrials = s1(1) / M;	
		for m = 1 : MM
			% % % shuffle A
			tmp1 = randperm(nTrials);
			tmp2 = bsxfun(@plus, tmp1 * M - M, (1 : M)');
			A1 = A(tmp2(:), :);
			dA1 = ~isnan(size(A));
			A1(dA1) = 0;
			
			% % % compute data again
			pj = zeros(s1(2), s2(2), N) + nan;
			for ib = 1 : N
				X = 
				sn = (da == ib) * dA1;
				sd = (da == ib) * A1;
				pj(:, :, ib) = sd ./ sn;
			end
			pm = bsxfun(@rdivide, pj, sum(pj, 3));
			pac_r(:, :, m) = 1 + sum(pm .* log(pm), 3) ./ log(N);
		end
	end
	
end %end of function

