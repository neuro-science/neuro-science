%% function #1
function [c, cp] = meg_cluster4cluster(d, ed, nb, pth, nth, nt)
	% % % written 08/05/2018 by wp : cluster for cluster
	
	% % % check inputs
	N = size(ed, 1);
	n = size(d, 3);
	if pth < 1
		pth = mat_tThreshFromP(n, pth);
	end
	
	% % % topology check
	[tp, ft] = tp_compute(ed, nb, N);
	fprintf('topo was processed @%04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
	
	% % % positive connection
	dd = d(:, 1, :) - d(:, 2, :);
	t = mat_tValue4Arrays(dd, 3);
	c1 = find(t > pth & ft);
	c2 = myFilt(c1, tp, ft, nth);
	[cc1, n1] = theCluster(c2, tp);
	fprintf('data+ was processed @%04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
	
	% % % negative connection
	c1 = find(-t > pth & ft);
	c2 = myFilt(c1, tp, ft, nth);
	[cc2, n2] = theCluster(c2, tp);
	fprintf('data- was processed @%04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
	
	% % % merge
	cc = [cc1; cc2];
	nc = n1 + n2;
	clear c1 c2 cc1 cc2 n1 n2;
	sz = zeros(nc, 1);
	for k = 1 : nc
		if ~isempty(cc{k})
			sz(k) = sum(t(cc{k}));
		end
	end
	[x1, x2] = sort(abs(sz), 'descend');
	idx = x2(logical(x1));
	cc = cc(idx);
	sz = sz(idx);
	fprintf('\nmax size : %f.\n\n', sz(1));
	
	% % % statistics
	s0 = [ones(n, 1); ones(n, 1) - 2];
	H0 = zeros(nt, 1);
	parfor it = 1 : nt
		s1 = randperm(n*2);
		s = s0(s1(1:n));
		dt = bsxfun(@mtimes, dd, permute(s, [2 3 1]));
		tt = mat_tValue4Arrays(dt, 3);
		c1 = find(tt > pth & ft);
		c2 = myFilt(c1, tp, ft, nth);
		[cc1, n1] = theCluster(c2, tp);
		if n1 > 0
			sz1 = zeros(n1, 1);
			for k = 1 : n1
				if ~isempty(cc1{k})
					sz1(k) = sum(tt(cc1{k}));
				end
			end
			H0(it) = max(abs(sz1));
		end
		fprintf('%.2f @ %03d\n', H0(it), it);
	end
	
	% % % sort data
	H = sort(H0, 'descend');
	th = H(floor(0.05 * nt));
	idx = abs(sz) > th;
	cc = cc(idx);
	c = cell(size(cc));
	for ic = 1 : length(cc)
		c{ic}.nd = cc{ic};
		c{ic}.sz = sz(ic);
	end
	cp.H = H;
	cp.th = th;
end % end of function

%% reciprocal 
function [y, fg] = theRcp(c, x, y, fg, tp)
	nn = numel(x);
	x2 = [];
	for k = 1 : nn
		if fg(x(k))
			fg(x(k)) = false;
			y = [y; x(k)];
			[x1, x2] = ismember(tp{c(x(k))}, c);
			x2 = x2(logical(x2));
			if ~isempty(x2)
				[y, fg] = theRcp(c, x2, y, fg, tp);
			end
		end
	end
end


%% cluster 
function [cc, nc] = theCluster(c, tp)
	nn = numel(c);
	fg = true(nn, 1);
	nc = 0;
	cc = cell(2, 1);
	for ie = 1 : nn
		y = [];
		[y, fg] = theRcp(c, ie, y, fg, tp);
		if ~isempty(y)
			nc = nc + 1;
			cc{nc} = sort(c(y));
		end
	end
end

%% filter
function c2 = myFilt(c1, tp, ft, nth)
	c2 = false(size(c1));
	for ik = 1 : length(c1)
		c2(ik) = sum(ismember(tp{c1(ik)}, c1)) ./ ft(c1(ik)) >= nth;
	end
	c2 = c1(c2);
end

%% the topo definition
function [tp, ft] = tp_compute(ed, nb, N)
	tp = cell(N, 1);
	ft = true(N, 1);
	for ie = 1 : N
		tmp1 = abs(bsxfun(@minus, ed(:, 1), ed(ie, 1)));
		tmp2 = abs(bsxfun(@minus, ed(:, 2), ed(ie, 2)));
		tmp3 = sum(abs(bsxfun(@minus, ed(:, 3:4), ed(ie, 3:4))), 2);
		tmp11 = find((~tmp2) & (~tmp3));
		s1 = tmp11(ismember(ed(tmp11, 1), nb{ed(ie, 1)}));
		tmp22 = find((~tmp1) & (~tmp3));
		s2 = tmp22(ismember(ed(tmp22, 1), nb{ed(ie, 2)}));
		s3 = find((~tmp1) & (~tmp2) & (tmp3 == 1));
		clear tmp*;
		tp{ie} = [s1; s2; s3];
		ft(ie) = numel(tp{ie});
	end
end

