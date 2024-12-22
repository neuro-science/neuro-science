function r = mat_mvpa4erf(x, y, np, na)
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
	
	for ip = np : -1 : 1
		fprintf('o');
		for ic = 2 : -1 : 1
			tmp = randperm(ntrls(ic));
			for it = n_(ic) : -1 : 1
				erf(n_(1) * (ic - 1) + it, :, :, ip) = mean(xx{ic}(tmp((it - 1) * na + 1 : it * na), :, :), 1);
			end
			for iv = sz(3) : -1 : 1
				s_(:, :, iv, 2, ip) = mat_mnn_cov(erf(n_(1) + 1 : n_(2), :, iv, ip));
				s_(:, :, iv, 1, ip) = mat_mnn_cov(erf(1 : n_(1), :, iv, ip));
			end
		end
	end
	fprintf('\n');
	clear xx;
	s = mean(mean(mean(s_, 3), 4), 5)^-0.5;
	if ~isreal(s)
		fprintf('complex s here!\n');
		s = real(s);
	end
	erf = reshape(permute(erf, [1 4 2 3]), [sum(n_)*np, sz(2), sz(3)]);
	for iv = sz(3) : -1 : 1
		erf(:, :, iv) = erf(:, :, iv) * s;
		fprintf('*');
	end
	fprintf('\n');
	clear s s_;
	
	% % % separate training and test trials		
	y = repmat([ones(n_(1), 1); 2 * ones(n_(2), 1)], [1 np]);
	y = y(:);
	ir = false(n_(1) + n_(2), np);
	ir([1 : ceil(n_(1)/2), n_(1) + 1 : n_(1) + ceil(n_(2)/2)], :) = true;
	ir = ir(:);
	is = ~ir;

	r = nan(sz(3), 1);
	for iv = 1 : sz(3)
% 	for iv = sz(3) : -1 : 1
		fprintf('.');
		model_svm = fitcsvm(erf(ir, :, iv), y(ir));
		r(iv) = mean(predict(model_svm, erf(is, :, iv)) == y(is));	
	end	
	fprintf('\n');
	clear erf model_svm;
	
end
