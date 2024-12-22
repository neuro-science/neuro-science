function [n, r] = mat_mvpa4vx(x, y, th, nr)
	% % % written by wp @22/03/2019 
	% % % mvpa for one voxel data
	% % % inputs:
	% % % x - data (trial x time)
	% % % y - label (trial x 1)
	% % % th - threshold of correct ratio
	% % % n1 - number of permutations
	% % % n2 - number of averaged trials at beginning
	% % % outputs:
	% % % n - minimum number of averages needed to be above threshold
	% % % th - thresholds obtained	
	%% 1. parse inputs	
	if nargin < 4 || isempty(n1)
		nr = 100;
	end

	if nargin < 3 || isempty(th)
		th = 0.85;
	end

	if nargin < 2
		fprintf('We need at least two inputs: data(trial X time) and label(trial X 1)!\n');
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
				xx{ic} = x(idx, :);
			end
		end
		clear tmp x y idx;
	end
	
	%% 3. work with stair-case like threshold finding
	mr = [0 0 0];
	nM = floor(min(ntrls)/2);
	na = 1;
	while mean(mr) < th && na <= nM 
		r = myTrainTest(xx, ntrls, nr, na);
		fprintf('current n = %3d r = %6.2f.\n', na, r);
		mr = [mr(2:3), r];
		na = na + 1;
	end
	n = na - 1;			
	r = mean(mr);
	clear xx;
end

function r = myTrainTest(x, n, nr, n2)
	
	if n2 == 1 && nr > 1
		fprintf('Only one trial for averaging, no permutations needed!\n');
		nr = 1;
	end
	
	rng(10);
	for ip = nr : -1 : 1
		% % % prepare data		
		for ic = 2 : -1 : 1
			n_(ic) = floor(n(ic) / n2);
		end
		for ic = 2 : -1 : 1
			tmp = randperm(n(ic));
			for it = n_(ic) : -1 : 1
				erf(n_(1) * (ic - 1) + it, :, ip) = mean(x{ic}(tmp((it - 1) * n2 + 1 : it * n2), :), 1);
			end
		end
		% % % noise process		
	end
	erf = reshape(permute(erf, [1 3 2]), sum(n_)*nr, []);
	s = (mat_mnn_cov(erf))^-0.5;
	if isreal(s)
		erf = erf * s;
	else
		erf = erf * real(s);
		fprintf('complex!\n');
	end
	
	
	% % % separate training and test trials		
	y = repmat([ones(n_(1), 1); 2 * ones(n_(2), 1)], [1 nr]);
	y = y(:);
	ir = false(n_(1) + n_(2), nr);
	ir([1 : ceil(n_(1)/2), n_(1) + 1 : n_(1) + ceil(n_(2)/2)], :) = true;
	ir = ir(:);
	is = ~ir;
	
	% % % do it		
	[size(y(ir)) size(y(is))]
	model_svm = fitcsvm(erf(ir, :), y(ir));
	r = mean(predict(model_svm, erf(is, :)) == y(is));	
	clear erf model_svm y x n*;
end






% % % % % % % % % % older versions % % % % % % % % % %
% function r = mat_mvpa4vx(x, y, np, na)
% % % written by wp @22/03/2019 
% % % mvpa for one voxel data
% % % inputs:
% % % x - data (trial x time x chan)
% % % y - label (trial x 1)
% % % np - number of permutations
% % % na - number of averaged trials
% % % outputs:
% % % n - minimum number of averages needed to be above threshold
% % % th - thresholds obtained	


% function [n, r] = mat_mvpa4vx(x, y, th, n1, n2)
% % % % written by wp @22/03/2019 
% % % % mvpa for one voxel data
% % % % inputs:
% % % % x - data (trial x time)
% % % % y - label (trial x 1)
% % % % th - threshold of correct ratio
% % % % n1 - number of permutations
% % % % n2 - number of averaged trials at beginning
% % % % outputs:
% % % % n - minimum number of averages needed to be above threshold
% % % % th - thresholds obtained	
% 
% 	%% 1. parse inputs
% 	if nargin < 5 || isempty(n2)
% 		n2 = 5;
% 	end
% 	
% 	if nargin < 4 || isempty(n1)
% 		n1 = 100;
% 	end
% 
% 	if nargin < 3 || isempty(th)
% 		th = 0.85;
% 	end
% 
% 	if nargin < 2
% 		fprintf('We need at least two inputs: data(trial X time) and label(trial X 1)!\n');
% 		return;
% 	else
% 		tmp = unique(y);
% 		if numel(unique(y)) ~= 2 
% 			fprintf('We currently support only two levels!\n');
% 			return;
% 		elseif numel(y) ~= size(x, 1)
% 			fprintf('We need consistent inputs: data(trial X time) and label(trial X 1)!\n');
% 			return;
% 		else
% 			for ic = 2 : -1 : 1
% 				idx{ic} = find(y == tmp(ic));
% 				ntrls(ic) = numel(idx{ic});
% 			end
% 			if min(ntrls) < n2*2
% 				fprintf('We need at least %d trials to perform the training!\n', n2*2);
% 				return;
% 			else
% 				for ic = 2 : -1 : 1
% 					xx{ic} = x(idx{ic}, :);
% 				end
% 			end
% 		end
% 		clear tmp x y idx;
% 	end
% 	
% 	
% 	%% 2. try with extreme cases first
% 	nM = floor(min(ntrls)/2);
% 	% % % 2.1 average	for half of the trials
% 	rM = myTrainTest(xx, ntrls, 10, nM);
% 	fprintf('Max = %6.2f.\n', rM);
% 	if rM < 0.51
% 		r = rM;
% 		n = nM;
% 		return;
% 	end
% 	% % % 2.2 without average	
% 	rm = myTrainTest(xx, ntrls, 1, 1);
% 	fprintf('Min = %6.2f.\n', rm);
% 	if rm > 0.99
% 		r = rm;
% 		n = 1;
% 		return;
% 	end
% 	
% 	%% 3. work with stair-case like threshold finding
% 	dr = [];
% 	while 1
% 		r = myTrainTest(xx, ntrls, n1, n2);
% 		fprintf('current n = %3d r = %6.2f.\n', n2, r);
% 		if isempty(dr)
% 			dr = r >= th;
% 			r_ = r;
% 		else
% 			dr_ = dr;
% 			dr = r >= th;
% 			if dr ~= dr_ 
% 				if dr
% 					break;
% 				else
% 					r = r_;
% 					n2 = n2 + 1;
% 					break;
% 				end
% 			end
% 		end
% 		n2 = n2 + 1 - 2 * dr;
% 		if n2 >= nM || n2 <= 1
% 			fprintf('cn = %3d r = %6.2f.\n', n2, r);
% 			break;
% 		end
% 	end
% 	n = n2;			
% 	clear xx;
% 	
% end
% 
% function r = myTrainTest(x, n, n1, n2)
% 	
% 	if n2 == 1 && n1 > 1
% 		fprintf('Only one trial for averaging, no permutations needed!\n');
% 		n1 = 1;
% 	end
% 	
% 	rng(10);
% 	for ip = n1 : -1 : 1
% 		% % % prepare data		
% 		for ic = 2 : -1 : 1
% 			n_(ic) = floor(n(ic) / n2);
% 		end
% 		for ic = 2 : -1 : 1
% 			tmp1 = randperm(n(ic));
% 			for it = n_(ic) : -1 : 1
% 				erf(n_(1) * (ic - 1) + it, :, ip) = mean(x{ic}((it - 1) * n2 + 1 : it * n2, :), 1);
% 			end
% 		end
% 		% % % noise process		
% 	end
% 	erf = reshape(permute(erf, [1 3 2]), sum(n_)*n1, []);
% 	s = (cov1para(erf))^-0.5;
% 	if isreal(s)
% 		erf = erf * s;
% 	else
% 		erf = erf * real(s);
% 		fprintf('complex!\n');
% 	end
% 	
% 	
% 	% % % separate training and test trials		
% 	y = repmat([ones(n_(1), 1); 2 * ones(n_(2), 1)], [1 n1]);
% 	y = y(:);
% 	ir = false(n_(1) + n_(2), n1);
% 	ir([1 : ceil(n_(1)/2), n_(1) + 1 : n_(1) + ceil(n_(2)/2)], :) = true;
% 	ir = ir(:);
% 	is = ~ir;
% 	
% 	% % % do it		
% 	model_svm = fitcsvm(erf(ir, :), y(ir));
% 	r = mean(predict(model_svm, erf(is, :)) == y(is));	
% 	clear erf model_svm y x n*;
% end
% 
% function r = myTrainTest(x, n, n1, n2)
% 	
% 	if n2 == 1 && n1 > 1
% 		fprintf('Only one trial for averaging, no permutations needed!\n');
% 		n1 = 1;
% 	end
% 	
% 	rng(10);
% 	for ip = n1 : -1 : 1
% 		% % % prepare data		
% 		for ic = 2 : -1 : 1
% 			n_(ic) = floor(n(ic) / n2);
% 		end
% 		for ic = 2 : -1 : 1
% 			tmp1 = randperm(n(ic));
% 			for it = n_(ic) : -1 : 1
% 				erf(n_(1) * (ic - 1) + it, :, ip) = mean(x{ic}((it - 1) * n2 + 1 : it * n2, :), 1);
% 			end
% 		end
% 		% % % noise process		
% 		s(:, :, 2, ip) = mat_mnn_cov(erf(n_(1) + 1 : n_(1) + n_(2), :, ip));
% 		s(:, :, 1, ip) = mat_mnn_cov(erf(1 : n_(1), :, ip));
% 	end
% 	s_ = mean(mean(s, 3), 4)^-0.5;
% 	if isreal(s_)
% 		erf = reshape(permute(erf, [1 3 2]), sum(n_)*n1, []) * s_;
% 	else
% 		s1 = sumsqr(real(s_));
% 		s2 = sumsqr(imag(s_));
% 		fprintf('complex s_ observed: %6.2f%%!\n', s2/s1);
% 		erf = reshape(permute(erf, [1 3 2]), sum(n_)*n1, []) * real(s_);
% 	end
% 	clear s s_;
% 	
% 	% % % separate training and test trials		
% 	y = repmat([ones(n_(1), 1); 2 * ones(n_(2), 1)], [1 n1]);
% 	y = y(:);
% 	ir = false(n_(1) + n_(2), n1);
% 	ir([1 : ceil(n_(1)/2), n_(1) + 1 : n_(1) + ceil(n_(2)/2)], :) = true;
% 	ir = ir(:);
% 	is = ~ir;
% 	
% 	% % % do it		
% 	model_svm = fitcsvm(erf(ir, :), y(ir));
% 	r = mean(predict(model_svm, erf(is, :)) == y(is));	
% 	clear erf model_svm y x n*;
% end