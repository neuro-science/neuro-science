function [r, r0] = mat_mvpa4avg(x, y, np, na, wFlag, tr, nr)
% % % updated by wp @28/03/2019 for permutation of null hypothesis
% % % updated by wp @26/03/2019 for improvement of noise normalization and flexible testing number
% % % written by wp @22/03/2019 
% % % mvpa for single subjects' data
% % % inputs:
% % % x - data (trial x time x chan)
% % % y - label (trial x 1)
% % % np - number of permutations
% % % na - number of averaged trials
% % % wFlag - whether whitening is needed
% % % tr - ratio/number of test trials per condition
% % % nr - number of permutation of null hypothesis for each permutation
% % % outputs:
% % % r - ratio of correct prediction
% % % r0 - ratio of correct prediction for null

	%% 1. parse inputs	
	if nargin < 6 || isempty(tr)
		rFlag = true;
		tr = 0.5;
	elseif tr > 0 && abs(mod(tr, 1)) < 1e-5
		rFlag = false;
	elseif tr > 0 && tr < 1
		rFlag = true;
	else
		fprintf('Please provide the 6th inputs as ratio (0~1) or number(positive integer) of test trials!\n');
		return;
	end

	if nargin < 5 || isempty(wFlag)
		wFlag = false;
	end

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
	
	%% 2. prepare id of labels etc
	rng(10);
	% % % prepare data		
	for ic = 2 : -1 : 1
		n_(ic) = floor(ntrls(ic) / na);
	end
	nn = sum(n_);
	
	% % % separate training and test trials		
	y = [ones(n_(1), 1); 2 * ones(n_(2), 1)];
	is = false(nn, 1);
	if rFlag
		xi1 = max(floor(n_(1)*tr), 1);
		xi2 = max(floor(n_(2)*tr), 1);
		is([1:xi1, n_(1)+1:n_(1)+xi2]) = true;
		clear xi1 xi2;
	else
		is([1:tr, n_(1)+1:n_(1)+tr]) = true;
	end
	ir = ~is;

	%% 3. do data preparison and testing for each permutation
	for ip = np : -1 : 1
		% % % prepare data		
		for ic = 2 : -1 : 1
			tmp = randperm(ntrls(ic));
			for it = n_(ic) : -1 : 1
				erf(n_(1) * (ic - 1) + it, :, :) = nanmean(xx{ic}(tmp((it - 1) * na + 1 : it * na), :, :), 1);
			end
		end
		% % % noise normalization by whitening
		if wFlag
			for iv = sz(3) : -1 : 1
				s_(:, :, iv, 2) = mat_mnn_cov(erf(n_(1) + 1 : nn, :, iv));
				s_(:, :, iv, 1) = mat_mnn_cov(erf(1 : n_(1), :, iv));
			end
			s = mean(mean(s_, 3), 4)^-0.5;
			if ~isreal(s)
				fprintf('x');
				s = real(s);
			end
			for iv = sz(3) : -1 : 1
				erf(:, :, iv) = erf(:, :, iv) * s;
			end
			clear s s_;
		end
		fprintf('-');
		% % % train and test        
		for iv = sz(3) : -1 : 1
			yy1 = y(ir);
			yy2 = y(is);
			ny1 = numel(yy1);
			x1 = erf(ir, :, iv);
			x2 = erf(is, :, iv);
			model_svm = fitcsvm(x1, yy1);
			rr(iv) = mean(predict(model_svm, x2) == yy2);	
			parfor ix = 1 : nr
				yy = yy1(randperm(ny1));
				model_svm = fitcsvm(x1, yy);
				r_(ix) = mean(predict(model_svm, x2) == yy2);	
			end
			r0_(:, iv) = r_;
			clear r_;
		end	
		r(:, ip) = rr;
		r0(:, :, ip) = r0_;
		clear erf rr r0_;
		fprintf('o');
	end
	clear xx;
	
end
