function [r, r0] = mat_mvpa4grp(x, wFlag, nr)
% % % written by wp @28/03/2019
% % % mvpa for group of averages
% % % inputs:
% % % x - data (sub x time x chan, 2), dimension time is used for decoding
% % % wFlag - whether whitening is needed
% % % nr - number of permutation of null hypothesis for each permutation
% % % outputs:
% % % r - ratio of correct prediction
% % % r0 - ratio of correct prediction for null

	%% 1. parse inputs	
	if nargin < 1
		fprintf('We need at least one input: data(sub X observeD X otherD X 2)!\n');
		return;
	else
		[nsb, nob, not, tmp] = size(x);
		if tmp ~= 2
			fprintf('We can only process two conditions for now!\n');
			return;
		end
		clear tmp;
	end
	
	if nargin < 2 || isempty(wFlag)
		wFlag = true;	%default with noise normalization
	end

	if nargin < 3 || isempty(nr)
		nr = 1000;
	end
	
	%% 2. prepare data, id of labels etc
	% % % labels	
	y = ones(nsb, 1) * [1 2];
	y = y(:);
	s1 = ~repmat(eye(nsb), [2, 1]);
	s2 = ~s1;
	
	% % % permutations estimation	
	rng(10);
	rn = nr / 2^(nsb-1);
	if rn > 1
		fprintf('The number of permutation is more than possible!\t');
	else
		fprintf('%6.2f%% of possible permutations!\t', rn * 100);
	end
	for ir = nr : -1 : 1
		y0(:, ir) = y(randperm(nsb*2));
	end
	% % % data noise normalization
	if wFlag
		for ic = 2 : -1 : 1
			for iv = not : -1 : 1
				s_(:, :, iv, ic) = mat_mnn_cov(x(:, :, iv, ic));
			end
		end
		ss = mean(mean(s_, 3), 4)^-0.5;
		if ~isreal(ss)
			fprintf('x\n');
			ss = real(ss);
		else
			fprintf('\n');
		end
		for ic = 2 : -1 : 1
			for iv = not : -1 : 1
				x1(:, :, iv, ic) = x(:, :, iv, ic) * ss;
			end
		end
		clear s_ ss;
	else
		x1 = x;
	end
	% % % reshape data	
	x = reshape(permute(x1, [1 4 2 3]), [nsb*2, nob, not]);	%[sb, time, vx, lb] - [sb lb time vx]
	clear x1;

	%% 3. do data fitting
	for is = nsb : -1 : 1
		y1 = y(s1(:, is));
		y2 = y(s2(:, is));
		y10 = y0(s1(:, is), :);
		y20 = y0(s2(:, is), :);
		for iv = not : -1 : 1
			x1 = x(s1(:, is), :, iv);
			x2 = x(s2(:, is), :, iv);
% 			theT = clock;fprintf('%03d@%02d:%02d ', not-iv+1, theT(4:5));clear theT;
			fprintf('.');
			model_svm = fitcsvm(x1, y1);
			r_(iv) = mean(predict(model_svm, x2) == y2);	
			rr = nan(nr, 1);
			parfor ir = 1 : nr
				y01 = y10(:, ir);
				y02 = y20(:, ir);
				model_svm = fitcsvm(x1, y01);
				rr(ir) = mean(predict(model_svm, x2) == y02);	
			end
			r0_(:, iv) = rr;
			clear rr x1 x2;
		end
		fprintf('\nDone for %d/%d @%04d-%02d-%02d %02d:%02d:%02d!\n', nsb-is+1, nsb, round(clock));
		r(:, is) = r_;
		r0(:, :, is) = r0_;
		fprintf('\n');
		clear r_ r0_;
	end
end
