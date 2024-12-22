function [c, th, H] = meg_cluster4ML(r, r0, para, th1, th2, n1, pth)
% % % written 03/05/2019 by wp : cluster for machine learning
% % % input r should be nvoxel x 1 or nt*nf x 1.
% % % input r0 should be N x nvoxel or N x nt*nf.
% % % outputs are clusters

	%% 1. check input
	if nargin < 7 || isempty(pth)
		pth = 0.05;
	end
	
	if nargin < 6 || isempty(n1)
		n1 = 3;
	end
	
	if nargin < 5 || isempty(th2)
		th2 = 0.5;
	end
	
	if nargin < 4 || isempty(th1)
		th1 = 0.5;
	end
	
	if nargin < 3
		fprintf('We need 3 inputs: data, null-distribution and neighborhood or tf size!\n');
		return;
	end
	
	if ~isvector(r)
		fprintf('The data should be a vector!\n');
		return;
	else
		n = numel(r);
	end
	
	sz = size(r0);
	if numel(sz) ~= 2
		fprintf('The null-distribution should be 2D!\n');
		return;
	elseif sz(2) ~= n && sz(1) ~= n
		fprintf('At least one dimension of the null-distribution should be same as data!\n');
		return;
	elseif sz(2) ~= n
		r0 = r0';
		sz = size(r0);
	end
	
	%% 2. clustering for data
	for iv = n : -1 : 1
% 		th3(iv, 1) = max(mat_distr2thresh (r0(:, iv), pth), th1);
		th3(iv, 1) = th1;
	end
	% % % binarize data	
	if iscell(para)	%neighborhood definition - spatial
		b = mySpatialFilter(r, th3, th2, para);
		cc = mySPCluster(b, para);
		H0 = zeros(sz(1), n1);
		parfor ik = 1 : sz(1)
			b = mySpatialFilter(r0(ik, :)', th3, th2, para);
			[tmp1, tmp2] = mySPCluster(b, para);
			tmp1 = [];
			% % % fill the results		
			szn = zeros(n1, 1);
			ssz = length(tmp2);
			if ssz >= n1
				szn(:, 1) = tmp2(1 : n1);
			elseif ssz > 0
				szn(1 : ssz, 1) = tmp2;
			end
			H0(ik, :) = szn;
			b = [];
			tmp2 = [];
		end
		[c, th, H] = sumClusterSP(H0, cc, n, pth);	
	elseif isvector(para) && numel(para) == 2	%tf definition - TF
		r = reshape(r, para);
		th3 = reshape(th3, para);
		r0 = reshape(r0, [para sz(1)]);
		b = mySpatialFilter(r, th3, th2);
		cc = myTFCluster(b);
		H0 = zeros(sz(2), n1);
		parfor ik = 1 : sz(1)
			b = mySpatialFilter(r0(:, :, ik), th3, th2);
			[tmp1, tmp2] = myTFCluster(b);
			tmp1 = [];
			% % % fill the results		
			szn = zeros(n1, 1);
			ssz = length(tmp2);
			if ssz >= n1
				szn(:, 1) = tmp2(1 : n1);
			elseif ssz > 0
				szn(1 : ssz, 1) = tmp2;
			end
			H0(ik, :) = szn;
			b = [];
			tmp2 = [];
		end
		[c, th, H] = sumClusterTF(H0, cc, para(1), para(2), pth);	
	end
	
end %end of function

%% 1. spatial filter
function b = mySpatialFilter(t, th1, th2, para)
	% % % thresh apply
	B = t > th1;
	if nargin < 4
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
	else
		for iv = numel(para) : -1 : 1
			B1(iv, 1) = sum(B(para{iv})) / numel(para{iv}) > th2;
		end
	end
	% % % output	
	b = B1 .* B .* t;
end

%% 2. TF cluster
function [c, sz] = myTFCluster(b)
	% % % get 1d vector of significant points
	[x1, x2, v] = find(b);
	ij = x1 > x2;
	x1(ij) = [];
	x2(ij) = [];
	v(ij) = [];

	% % % do it 	
	nds = length(x1);
	flag = true(nds, 1);
	ids = [x1, x2, v];
	nc = 0;
	for ig = 1 : nds
		if flag(ig)
			cid = [];
			[cid, flag] = myTFSearch(ids(:, [1 2]), ig, cid, flag);
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

%% 3. TF search
function [y, flag] = myTFSearch(ids, x, y, flag)
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
		[y, flag] = myTFSearch(ids, x2, y, flag);
	end
end

%% 4. spatial cluster
function [c, sz] = mySPCluster(b, para)
	% % % get 1d vector of significant points
	ids = find(b);
	if isempty(ids)
		c = [];
		sz = [];
		return;
	else
		flag = true(size(ids));
		nc = 0;
		for iv = 1 : numel(ids)
			cid = [];
			[cid, flag] = mySPSearch(ids, iv, cid, flag, para);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = numel(cid);
				c{nc}.ed(:, 1) = ids(cid);
				c{nc}.ed(:, 2) = b(ids(cid));
				c{nc}.sz = sum(c{nc}.ed(:, 2));
				sz(nc) = c{nc}.sz;
			end
		end
	end
	clear b flag ids iv cid para;
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

%% 5. spatial search
function [y, flag] = mySPSearch(ids, x, y, flag, para)
	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			x1 = para{ids(x(k))};
			[tmp1, tmp2] = ismember(x1, ids);
			x2 = [x2; tmp2(tmp1)];
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = mySPSearch(ids, x2, y, flag, para);
	end
end

%% 6. summarySP
function	[c, th, H] = sumClusterSP(H0, cc, n, pth)	
	% % % get the thresh
	[th, H] = mat_distr2thresh(H0, pth);

	% % % select the clusters
	c = cell(2, 1);
	ct = 0;
	for ik = 1 : numel(cc)
		if cc{ik}.sz > th
			ct = ct + 1;
			c{ct} = cc{ik};
			c{ct}.p0 = cluster_pValue(H, cc{ik}.sz);
			c{ct}.spt = zeros(n, 1);
			c{ct}.spt(c{ct}.ed(:, 1)) = c{ct}.ed(:, 2);
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


%% 7. summaryTF
function	[c, th, H] = sumClusterTF(H0, cc, nf, nt, pth)	
	% % % get the thresh
	[th, H] = mat_distr2thresh(H0, pth);

	% % % select the clusters
	c = cell(2, 1);
	ct = 0;
	for ik = 1 : numel(cc)
		if cc{ik}.sz > th
			ct = ct + 1;
			c{ct} = cc{ik};
			c{ct}.p0 = cluster_pValue(H, cc{ik}.sz);
			i = c{ct}.ed(:, 1);
			j = c{ct}.ed(:, 2);
			v = c{ct}.ed(:, 3);
			c{ct}.tft = full(sparse(i, j, v, nf, nt));
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

%% 8. pvalue
function p = cluster_pValue(H, sz)
	[y, I] = min(abs(H - sz));
	p = I ./ length(H);
	clear y I H sz;
end %end of function


