function c = meg_cluster1voxel(d, th1, th2, n1, n2, pth, fname)
% % % written 17/08/2018 by wp : cluster for one voxel
% % % input d should be in (f, t/f2, nsb) size

	% % % read size and initialize
	[nf, nt, ns] = size(d);
	% % % check threshold
	if th1 < 0.5 && th1 > 0
		th1 = mat_tThreshFromP(ns, th1);
	end
	% % % t-value	
	t = mat_tValue4Arrays(d);
	% % % find cluster	
	b = mySpatialFilter(-t, th1, th2);
	cc{2} = myFreeCluster(b);
	b = mySpatialFilter(t, th1, th2);
	cc{1} = myFreeCluster(b);
	save(fname, 'cc');

	%% initialize
	H0 = zeros(n1, n2);
	if n2 > 2^ns
		error('The number of subjects is too small!');
	elseif ns <= 20
		s0 = uint32(randperm(2.^ns));
		s1 = ismember(dec2bin(s0(1 : n2) - 1, ns), '1');
		clear s0;
	else
		s1 = rand(n2, ns) > 0.5;
	end
	
	%% get statistics
	parfor ts = 1 : n2
% 		for ts = 1 : n2
		tic;

		% % % prepare data		
		cz = bsxfun(@times, d, 2 * (permute(s1(ts, :), [1 3 2]) - 0.5));
		con = mat_tValue4Arrays(cz, 3, 0);
		b = mySpatialFilter(con, th1, th2);

		% % % do cluster
		[tmp1, sz0] = myFreeCluster (b);
		tmp1 = [];
		tmp3 = [];

		% % % fill the results		
		szn = zeros(n1, 1);
		ssz = length(sz0);
		if ssz >= n1
			szn(:, 1) = sz0(1 : n1);
		elseif ssz > 0
			szn(1 : ssz, 1) = sz0;
		end
		H0(:, ts) = szn;
		cz = [];
		con = [];
		b = [];
		% % % echo message		
		fprintf('Test #%03d done after %7.2f seconds.\n', ts, toc);
	end
	save(fname, 'H0', '-append');
	
	%% summary
	[c, th, H] = sumCluster(H0, cc, nf, nt, pth);	
	save(fname, 'c', 'th', 'H', '-append');
	
end %end of function

%% 1. spatial filter
function b = mySpatialFilter(t, th1, th2)
	% % % thresh apply
	B = t > th1;
	% % % TF neighbor check
	B4 = zeros([size(t), 4]);
	% % % shift array to get significant neighbors	
	B4(2 : end, :, 1) = B(1 : end - 1, : , :);
	B4(1 : end - 1, :, 2) = B(2 : end, : , :);
	B4(:, 2 : end, 3) = B(:, 1 : end - 1, :);
	B4(:, 1 : end - 1, 4) = B(:, 2 : end, :);
	% % % possible neighbors	
	B0 = zeros(size(t)) + 4;
	B0(1, :) = B0(1, :) - 1;
	B0(:, 1) = B0(:, 1) - 1;
	B0(end, :) = B0(end, :) - 1;
	B0(:, end) = B0(:, end) - 1;
	% % % apply thresh	
	B1 = sum(B4, 3)./B0 >= th2;
	% % % output	
	b = B1 .* B .* t;
end

%% 2. free cluster
function [c, sz] = myFreeCluster(b)
	% % % get 1d vector of significant points
	[x1, x2, v] = find(b);
% 	plot(x1, x2, '.');
	ij = x1 > x2;
	x1(ij) = [];
	x2(ij) = [];
	v(ij) = [];
% 	hold on;
% 	plot(x1, x2, 'o');

	% % % do it 	
	nds = length(x1);
	flag = true(nds, 1);
	ids = [x1, x2, v];
	nc = 0;
	for ig = 1 : nds
		if flag(ig)
			cid = [];
			[cid, flag] = myFreeSearch(ids(:, [1 2]), ig, cid, flag);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :);
				c{nc}.sz = sum(c{nc}.ed(:, 3));
				sz(nc) = c{nc}.sz;
			end
		end
	end
	clear x1 x2 v b nds flag ids ig cid;
	if nc
		[y, I] = sort(abs(sz), 'descend');
		c = c(I);
		sz = sz(I);
		clear I y;
	else
		c = [];
		sz = [];
	end
end

%% 3. free search
function [y, flag] = myFreeSearch(ids, x, y, flag)
	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			s1 = sum(bsxfun(@eq, ids, ids(x(k), :)), 2) == 1;	%1dim same
			s2 = sum(bsxfun(@eq, ids, ids(x(k), :)-1), 2) == 1;	%1dim -1
			s3 = sum(bsxfun(@eq, ids, ids(x(k), :)+1), 2) == 1;	%1dim +1
			s = find(s1 & (s2 | s3));
			x2 = [x2; s];
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = myFreeSearch(ids, x2, y, flag);
	end
end

%% 4. summary
function	[c, th, H] = sumCluster(H0, cc, nf, nt, pth)	
	% % % get the thresh
	[th, H] = mat_distr2thresh(H0, pth);
	

	% % % select the clusters
	c = cell(2, 1);
	ct = 0;
	for k = 1 : 2
		for ik = 1 : numel(cc{k})
			if cc{k}{ik}.sz > th
				ct = ct + 1;
				c{ct} = cc{k}{ik};
				c{ct}.p0 = cluster_pValue(H, cc{k}{ik}.sz);
				c{ct}.sgn = 2 * (1.5 - k);
				i = c{ct}.ed(:, 1);
				j = c{ct}.ed(:, 2);
				v = c{ct}.ed(:, 3);
				c{ct}.tft = full(sparse(i, j, v, nf, nt));
			end
		end
	end
	
	% % %  sort the clusters
	if ct
		for k = ct : -1 : 1
			sz(k) = c{k}.sz;
		end
		[tmp, I] = sort(abs(sz), 'descend'); %modeified 09/05/2017
		c = c(I);
		clear tmp I sz;
	else
		c = [];
	end
	clear ct cc k ik i j v;
	fprintf('=========Data sorted on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));
end

%% 5. 
function p = cluster_pValue(H, sz)
	[y, I] = min(abs(H - sz));
	p = I ./ length(H);
	clear y I H sz;
end %end of function

%% unknown previous version, a bit odd%% 
% % % % function [tn, ts] = meg_cluster1voxel(t, nb, th1, th2, N)
% % % % % % % written 17/08/2018 by wp : cluster for one voxel
% % % % % % % input t should be in (f, t, v) size
% % % % 
% % % % 	% % % read size and initialize
% % % % 	[nf, nt, nv] = size(t);
% % % % 	tn = zeros(nv, 1);
% % % % 	ts = zeros(size(t));
% % % % 	% % % check threshold
% % % % 	if nargin > 4 && ~isempty(N) && th1 < 0.5
% % % % 		th1 = mat_tThreshFromP(th1, N);
% % % % 	end
% % % % 	B = abs(t) > th1;
% % % % 	% % % TF neighbor
% % % % 	B4(:, :, :, 1) = B([1 1:end-1], :, :);
% % % % 	B4(:, :, :, 2) = B([2:end, end], :, :);
% % % % 	B4(:, :, :, 3) = B(:, [1 1:end-1], :);
% % % % 	B4(:, :, :, 4) = B(:, [2:end, end], :);
% % % % 	B4 = mean(B4, 4) >= 0.5;
% % % % 	for iv = 1 : nv
% % % % 		% % % spatial neighbor
% % % % 		B1 = mean(B(:, :, nb{iv}), 3) > th2;
% % % % 		B2 = B1 & B4(:, :, iv) & B1;
% % % % 		t1 = t(:, :, iv) .* B2;
% % % % 		tn(iv) = sum(abs(t1(:)));
% % % % 		ts(:, :, iv) = t1;
% % % % 		clear B1 B2 t1;
% % % % 	end
% % % % 
% % % % end %end of function

