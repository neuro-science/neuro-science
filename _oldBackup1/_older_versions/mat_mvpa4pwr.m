function r = mat_mvpa4pwr(x, y, np, na)
% % % written by wp @22/03/2019 
% % % mvpa for one voxel data
% % % inputs:
% % % x - data (trial x time x chan)
% % % y - label (trial x 1)
% % % np - number of permutations
% % % na - number of averaged trials
% % % outputs:
% % % r - ratio of correct prediction

	%% 1. parse inputs	
	if nargin < 4 || isempty(na)
		na = 20;
	end

	if nargin < 3 || isempty(np)
		np = 100;
	end

	if nargin < 2
		fprintf('We need at least two inputs: data(trial X time X chans) and label(trial X 1)!\n');
		return;
	else
		tmp = unique(y);
		if numel(unique(y)) ~= 2 
			fprintf('We currently support only two levels!\n');
			return;
		elseif numel(y) ~= size(x, 1)
			fprintf('We need consistent inputs: data(trial X time) and label(trial X 1)!\n');
			return;
		else
			for ic = 2 : -1 : 1
				idx = find(y == tmp(ic));
				ntrls(ic) = numel(idx);
				xx{ic} = x(idx, :, :);
			end
			sz = size(x);
		end
		clear tmp x y idx;
	end
	
	%% 2. work with stair-case like threshold finding
	rng(10);
	% % % prepare data		
	for ic = 2 : -1 : 1
		n_(ic) = floor(ntrls(ic) / na);
	end
	
	% % % separate training and test trials		
	y = [ones(n_(1), 1); 2 * ones(n_(2), 1)];
% 	ir = false(n_(1) + n_(2), 1);
% 	ir([1 : ceil(n_(1)/2), n_(1) + 1 : n_(1) + ceil(n_(2)/2)]) = true;
% 	is = ~ir;
	is = false(n_(1) + n_(2), 1);
	is([1, n_(1) + 1]) = true;
	ir = ~is;
	
	for ip = np : -1 : 1
% 		tic;
		% % % prepare data		
		fprintf('o ');
		for ic = 2 : -1 : 1
			tmp = randperm(ntrls(ic));
			for it = n_(ic) : -1 : 1
				erf(n_(1) * (ic - 1) + it, :, :) = nanmean(xx{ic}(tmp((it - 1) * na + 1 : it * na), :, :), 1);
			end
		end
% 		toc;tic;
		for iv = sz(3) : -1 : 1
			s_(:, :, iv, 2) = mat_mnn_cov(erf(n_(1) + 1 : sum(n_), :, iv));
			s_(:, :, iv, 1) = mat_mnn_cov(erf(1 : n_(1), :, iv));
		end
		s = mean(mean(s_, 3), 4)^-0.5;
		if ~isreal(s)
			s1 = sum(real(s(:)).^2);
			s2 = sum(imag(s(:)).^2);
			fprintf('complex s %7.2f%% here!\t', 100 * s1/s2);
			s = real(s);
		end
		for iv = sz(3) : -1 : 1
			erf(:, :, iv) = erf(:, :, iv) * s;
		end
% 		toc;tic;
		clear rr;
		for iv = sz(3) : -1 : 1
			model_svm = fitcsvm(erf(ir, :, iv), y(ir));
			rr(iv) = mean(predict(model_svm, erf(is, :, iv)) == y(is));	
		end	
		r(:, ip) = rr;
		clear s s_ erf;
% 		toc;
		fprintf('Done @%04d-%02d-%02d %02d:%02d:%02d!\n', round(clock));
	end
	clear xx;
	
end
