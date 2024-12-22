function [pac, pac_r, out] = meg_crossfreq_MI(a, A, N, M, MM, outRFlag)%, nPools
% % % 16/01/17	modified by wp, increase speed at price of memory 
% % % 16/01/17	modified by wp, increase speed at price of memory 
% % % 12/01/17	written by wp 
% % % compute modulation Index based on entropy measures
% % % 	a: angle of low frequency data (points x channels)
% % % 	A: Amplitude of high frequency data (points x channels)
% % %		N: number of bins within a cycle, default 12

	%% 1. check inputs
% 	if nargin < 7 || isempty(nPools)
% 		nPools = 4;
% 	end
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
	da = zeros(s1) + nan;
	for c1 = 1 : s1(2)
		[tmp1, da(:, c1)] = histc(a(:, c1), borders);
	end
	da = da';
	
	% % % bigger matrix
	X = false([N, size(da)]); % [bins, chan, data]
	for ib = 1 : N
		X(ib, :, :) = da == ib;
	end
	X = reshape(X, [N * s1(2), s1(1)]);
	
	% % % take the amp channel in
	dA = ~isnan(A);
	A(~dA) = 0;
	sn = X * double(dA);
	sd = X * double(A);
	sp = reshape(sd ./ sn, [N, s1(2), s2(2)]);
	pj = bsxfun(@rdivide, sp, sum(sp, 1));
	if nargin > 5 
		if outRFlag
			out = cell(2, 1);
			dPJ = pj;
			rPJ = zeros([N, s1(2), s2(2), MM]);
		else
			out = pj;
		end
	else
		outRFlag = false;
		out = [];
	end
	% % % results
	pac = squeeze(1 + sum(pj .* log(pj), 1) ./ log(N));
	clear sp sn sd dA borders da ib;

	
	%% 3. statistics 
	if randFlag % with permutation
		% % % initialize
		pac_r = zeros([s1(2), s2(2), MM]) + nan;
		nTrials = floor(s1(1) / M);	
		if MM > 30
% 			try
% 				thePool = parpool(nPools); 
% 				fprintf('\nIt seemed the parallel is OK!\n');
% 			catch ME
% 				fprintf('\nIt seemed the parallel is not working!\n');
% 				fprintf(ME.message);
% 				return;
% 			end
			parfor m = 1 : MM
				% % % shuffle A
				tmp1 = randperm(nTrials);
				tmp2 = bsxfun(@plus, tmp1 * M - M, (1 : M)');
				A1 = A(tmp2(:), :);
				dA1 = ~isnan(A1);
				A1(~dA1) = 0;

				% % % compute data again
				sn = X * double(dA1);
				sd = X * A1;
				sp = reshape(sd ./ sn, [N, s1(2), s2(2)]);
				pj = bsxfun(@rdivide, sp, sum(sp, 1));
				pac_r(:, :, m) = 1 + sum(pj .* log(pj), 1) ./ log(N);
				if outRFlag
					rPJ(:, :, :, m) = pj;
				end
				sp = []; sn = []; sd = []; pj = [];
			end
% 			delete(thePool);
			if outRFlag
				out{1} = dPJ;
				out{2} = rPJ;
				clear dPJ rPJ;
			end
		else
			for m = 1 : MM
				% % % shuffle A
				tmp1 = randperm(nTrials);
				tmp2 = bsxfun(@plus, tmp1 * M - M, (1 : M)');
				A1 = A(tmp2(:), :);
				dA1 = ~isnan(A1);
				A1(~dA1) = 0;

				% % % compute data again
				sn = X * double(dA1);
				sd = X * A1;
				sp = reshape(sd ./ sn, [N, s1(2), s2(2)]);
				pj = bsxfun(@rdivide, sp, sum(sp, 1));
				pac_r(:, :, m) = 1 + sum(pj .* log(pj), 1) ./ log(N);
				if outRFlag
					rPJ(:, :, :, m) = pj;
				end
				sp = []; sn = []; sd = []; pj = [];
			end
			if outRFlag
				out{1} = dPJ;
				out{2} = rPJ;
				clear dPJ rPJ;
			end
		end
	end
	
end %end of function

