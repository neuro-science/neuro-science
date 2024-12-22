function [c, r, pp] = meg_cluster4pac (conArray, neighborDefinition, conThresh, spfThresh, ...
	method, shiftPredictor, nClustersPerTest, nTestPerCluster, fname, nPools, clusterP, sumFlag)


% cluster4pac will find clusters based on the phase amplitude coupling array, 
% which defines connection strength between locations in various time. 
% Use as
%   [clusters, sizes] = cluster4pac (conArray, neighborDefinition, conThresh, spfThresh, method, shiftPredictor)
% 
% Input:
%		conArray				- a PhaseLoc x AmpLoc x Time array x numSubs x sti/control, defining the 
%									connection strength between two locations at certain time 
%		neighborDefinition	- a numGrid x 1 cell array, k-th element contains neigboring 
%									grids of grid k, including itself.
%		conThresh				- threshhold of connection strength, value in conArray 
%									above it will indicate a connection or else not. 
%		spfThresh				- spatial threshold filter will be applied to remove spurious connection.
%		method               - method of clustering, can be 'node' or 'edge' 
%		shiftPredictor			- same size as conArray, provide controls of	shift predictor
%		fname						- file name to store the results
% 
% Output:
%		clusters - a cell array, each element contains a cluster
%		sizes - accumulated values of strength
% 
% 
% Written by Peng Wang, Institute of Neurophysiology, UKE.
%
%   Modification History:
%   17/09/2017: first draft, adopted from the previous function of coupling cluster


	
	%% 1. parameters preparision
	fprintf('\n=========Cluster for PAC begins @%04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));
	% % % input paras and default
	if nargin < 2
		error('At least two inputs needed: connection array and neigbor definition');
	end
	
	if nargin < 3 || isempty(conThresh)
		pp.conThresh = 0.05;
	else
		pp.conThresh = conThresh;
	end
	
	if nargin < 4 || isempty(spfThresh)
		spfThresh = 0.5;	% 
	else
		pp.spfThresh = spfThresh;
	end
	
	if nargin < 5 || isempty(method)
		pp.method = 'node';
	else
		pp.method = method;
	end
	
	if nargin < 6 || isempty(shiftPredictor)
		pp.shiftFlag = false;
	else
		pp.shiftFlag = true;
	end
	
	if nargin < 7 || isempty(nClustersPerTest)
		pp.nClustersPerTest = 3;
	else
		pp.nClustersPerTest = nClustersPerTest;
	end
	
	if nargin < 8 || isempty(nTestPerCluster)
		pp.nTestPerCluster = 10;
	else
		pp.nTestPerCluster = nTestPerCluster;
	end
	
	if nargin < 9 || isempty(fname)
		pp.fname = '~/tmp4clusterRunning.mat';
	else
		pp.fname = fname;	
	end
	
	if nargin < 10 || isempty(nPools)
		pp.nPools = 4;
	else
		pp.nPools = nPools;
	end
		
	if nargin < 11 || isempty(clusterP)
		pp.clusterP = 0.05;
	else
		pp.clusterP = clusterP;
	end
		
	if nargin < 12 || isempty(sumFlag)
		pp.sumFlag = false;
	else
		pp.sumFlag = sumFlag;
	end
	
	
	% % % compute the paras and some check
	[pp.n, n1, pp.nTs, nSubs, n2] = size(conArray);
	
	tCon = mat_tThreshFromP(nSubs, pp.conThresh);
	t05 = mat_tThreshFromP(nSubs, 0.05, [], [], 0);

	if abs(pp.n - n1) > 0.1
		fprintf('\n Raw Data size error!\n');
		return;
	end
	if pp.shiftFlag 
		if max(abs(size(conArray) - size(shiftPredictor))) > 0.1
			fprintf('\n Shift Data size error!\n');
			return;
		end
		tShift = mat_tValue4Arrays(conArray(:, :, :, :, 1) - shiftPredictor(:, :, :, :, 1), 4);
		fShift1 = tShift > t05;
		fShift2 = fShift1;
	else
		fShift1 = true(size(conArray(:, :, :, 1)));
		fShift2 = fShift1;
	end
			
	if n2 == 2
		if pp.shiftFlag
			tShift = mat_tValue4Arrays(conArray(:, :, :, :, 2) - shiftPredictor(:, :, :, :, 2), 4);
			fShift2 = tShift > t05;
			shiftPredictor = shiftPredictor(:, :, :, :, 1) - shiftPredictor(:, :, :, :, 2);
			conArray = conArray(:, :, :, :, 1) - conArray(:, :, :, :, 2);
			tShift = mat_tValue4Arrays(abs(conArray) - abs(shiftPredictor), 4);
			fShift1 = fShift1 & (tShift > t05);
			fShift2 = fShift2 & (tShift > t05);
		else
			conArray = conArray(:, :, :, :, 1) - conArray(:, :, :, :, 2);
		end
	elseif n2 ~= 1
		fprintf('\n Data size error!\n');
		return;
	end
	
	% % % shortcut to run only summary
	if pp.sumFlag && exist(fname, 'file')
		cmd = input('\nThe output file was there already, what should we do?\n(y)es to load it and (n)o to ignore it:\n', 's');
		if strcmpi (cmd(1), 'y' )
			v = load(fname);
			if isfield(v, 'r') && isfield(v.r, 'H') && ~isempty(v.r.H)
				r = v.r; clear v;
				t = mat_tValue4Arrays(conArray(:, :, :, :, 1) - conArray(:, :, :, :, 2), 4);
				c = sum_clusters(r, pp, t);
				save(fname, 'c', 'pp', '-append');
				fprintf('\n We summarized the data and resaved!\n');
				return;
			end
		end
	end
	
	% % % test the file written and pools
	try
		save(fname, 'pp');
		fprintf('\nThe output file seems OK to write, checkup of parameters done!\n');
	catch ME
		fprintf('Output failed!');
		fprintf(ME.message);
		return;
	end
	
	%% 3. cluster +/-
	fprintf('\n%s =========Data computed in on %04d-%02d-%02d %02d:%02d:%02d=========\n', method, round(clock));
	t = mat_tValue4Arrays(conArray, 4, 0);
	[r.c{1}, r.nc1, r.sz1] = insFindCluster(t, tCon, neighborDefinition, pp.spfThresh, pp.method, fShift1);
	[r.c{2}, r.nc2, r.sz2] = insFindCluster(-t, tCon, neighborDefinition, pp.spfThresh, pp.method, fShift2);
	save(fname, 'r', 'pp', 't');
	fprintf('\n%s =========Data computed out on %04d-%02d-%02d %02d:%02d:%02d=========\n', method, round(clock));
	if ~isempty(r.sz1) && ~isempty(r.sz2)
		fprintf('\nMax size :%f \n', max(r.sz1(1), r.sz2(1)));
	elseif ~isempty(r.sz1)
		fprintf('\nMax size :%f \n', r.sz1(1));
	elseif ~isempty(r.sz2)
		fprintf('\nMax size :%f \n', r.sz2(1));
	else
		fprintf('\nEmpty size\n');
		return;
	end
	
	%% 4. statistical test
	% % % initialize
	H0 = zeros(pp.nClustersPerTest, pp.nTestPerCluster);
	if nSubs < 12.5
		s0 = randperm(2.^nSubs) - 1; %to avoid the nSubs+1 bit
		s1 = ismember(dec2bin(s0(1 : nTestPerCluster), nSubs), '1');
	else
		s1 = rand(nTestPerCluster, nSubs) > 0.5;
	end
	clear s0;
	
	if pp.nPools > 1
		% % % try the pools
		try
			thePool = parpool(pp.nPools); 
			fprintf('\nIt seemed the parallel is OK!\n');
		catch ME
			fprintf('\nIt seemed the parallel is not working!\n');
			fprintf(ME.message);
			return;
		end
		
		% % % do with parallel computing
		parfor ts = 1 : pp.nTestPerCluster
			tic;
			% % % prepare data	
			cz = bsxfun(@times, conArray, 2 * (permute(s1(ts, :), [1 3 4 2]) - 0.5));
			t1 = mat_tValue4Arrays(cz, 4, 0);

			% % % do clustering
			[tmp1, tmp2, sz0] = insFindCluster(t1, tCon, neighborDefinition, pp.spfThresh, pp.method, fShift1 & fShift2);
			tmp1 = [];
			tmp2 = [];
			cz = [];
			t1 = [];
			
			% % % fill the results		
			szn = zeros(pp.nClustersPerTest, 1);
			ssz = numel(sz0);
			if ssz >= pp.nClustersPerTest
				szn(:, 1) = sz0(1 : pp.nClustersPerTest);
				fprintf('Hx = %10.1f.\tTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
			elseif ssz > 0
				szn(1 : ssz, 1) = sz0;
				fprintf('Hx = %10.1f.\tTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
			else
				fprintf('Hx = \tnan.\tTest #%03d done after %7.2f seconds.\n', ts, toc);
			end
			H0(:, ts) = szn;

% 			% % % echo message
% 			try 
% 				fprintf('Hx = %10.1f.\tTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
% 			catch ME
% 				fprintf([ME.message, 'Test #%03d done after %7.2f seconds.\n'], ts, toc);
% 			end			
		end
		delete(thePool);
	else
		% % % do with parallel computing
		for ts = 1 : pp.nTestPerCluster
			tic;
			% % % prepare data	
			cz = bsxfun(@times, conArray, 2 * (permute(s1(ts, :), [1 3 4 2]) - 0.5));
			t1 = mat_tValue4Arrays(cz, 4, 0);
			cz = []; 
			toc;tic;
			% % % do clustering
			[tmp1, tmp2, sz0] = insFindCluster(t1, tCon, neighborDefinition, pp.spfThresh, pp.method, fShift1 & fShift2);
			tmp1 = [];
			tmp2 = [];
			t1 = [];

			% % % fill the results		
			szn = zeros(pp.nClustersPerTest, 1);
			ssz = length(sz0);
			if ssz >= pp.nClustersPerTest
				szn(:, 1) = sz0(1 : pp.nClustersPerTest);
				fprintf('Hx = %10.1f.\tTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
			elseif ssz > 0
				szn(1 : ssz, 1) = sz0;
				fprintf('Hx = %10.1f.\tTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
			else
				fprintf('Hx = \tnan.\tTest #%03d done after %7.2f seconds.\n', ts, toc);
			end
			H0(:, ts) = szn;
			
% 			% % % echo message
% 			try 
% 				fprintf('Hx = %10.1f.\tTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
% 			catch ME
% 				fprintf([ME.message, '\tTest #%03d done after %7.2f seconds.\n'], ts, toc);
% 			end			
		end
	end
	
	% % % get the threshold
	r.H0 = H0'; %[], e.g. [3, 10] -> [10, 3];
	[r.th, r.H] = mat_distr2thresh(r.H0, pp.clusterP);

	% % % save the data
	save(fname, 'r', '-append');
	
	%% 5. summarize the results
	
	% % % select the clusters
	c = sum_clusters(r, pp, t);
	
	% % % save results
	save(fname, 'c', '-append');
	
	% % % cleanup and notify
	fprintf('=========Data sorted and saved on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));
	
end



%% subfunctions 1: find clusters 
function [c, nc, sz] = insFindCluster(t, t0, nb, th, method, ctrl)

	% % % binarize
	bA = t > t0;

	% % % control?
	if nargin > 5
		bA = bA & ctrl;
	end
	
	% % % spatial filter
	bA = spatialFilter4Clusters(bA, nb, th);

	% % % find cluster
	switch lower(method(1))
		case 'n'
			[c, nc] = cluster4Node(bA);
		case 'e'
			[c, nc] = cluster4Edge(bA, nb);
		otherwise
			error('Current supported methods are <node> and <edge> only!');
	end
	
	% % % get the size of cluster and sort 
	if nc > 0
		tmp = zeros(nc, 1);
		nn = zeros(nc, 1);
		for k = 1 : nc
			tmp(k) = sum(t(c{k}.id));
			c{k}.sz = tmp(k);
		end
		[sz, I] = sort(tmp, 'descend');
		c = c(I);
	else
		sz = [];
	end
	
	% % % clean up
	clear t bA nb;
end

%% subfunctions 2: spatial filters
function out = spatialFilter4Clusters(in, nb, th)
	
	% % % paras setting
	sz = size(in);	%expect(node, node, time)
	if length(sz) < 3
		sz(3) = 1;
	end
	n = length(nb);
	
	% % % possible neighbors for each edge
	for k = n : -1 : 1
		np1d(k) = numel(nb{k});
	end
% 	np_edge = bsxfun(@plus, np1d, np1d') - 2;
	np_edge = bsxfun(@plus, np1d, np1d');
% % % 	np_edge = zeros(n); %[n, n, 1]
% % % 	for k1 = 1 : n
% % % 		for k2 = 1 : n
% % % 			if k1 <= k2
% % % 				np_edge(k1, k2) = numel(nb{k1}) + numel(nb{k1}) - 2;
% % % 			else
% % % 				np_edge(k1, k2) = np_edge(k2, k1);
% % % 			end
% % % 		end
% % % 	end
	
	% % % possible neighbors for each time
	np_time = zeros(sz(3), 1) + 2;
	np_time(1) = 1;
	np_time(end) = 1;
	np_time = permute(np_time, [2 3 1]); %[T, 1, 1] - [1 1 T]
	
	% % % real number of Neighbors
	if sz(3) == 1
		nTime = ones(size(in));
	else
		nTime = timeFilter(in);
	end
	nSpace = spaceFilter(in, nb);
	
	% % % sum together and apply filter	
	np_All = bsxfun(@plus, np_edge, np_time); %[n, n, T]
	out = in & ((nTime + nSpace) > th * np_All); %[n n T]
end

%% subfunction 3 node clusters
function [c, nc] = cluster4Node(in)
	% % % input check
	sz =size(in); %[n n T]
	[x1, x2, x3] = ind2sub(sz, find(in));
	nds = length(x1);
	fg = true(nds, 1);
	ids = [x1, x2, x3]; %[n n T]
	nc = 0;
	c = [];
	for ig = 1 : nds
		if fg(ig)
			cid = [];
			[cid, fg] = searchClusterLoopSparseNode(ids, ig, cid, fg);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :); % [m x 3(n n T)]
				c{nc}.nd1 = unique(ids(cid, 1));
				c{nc}.nn1 = length(c{nc}.nd1);
				c{nc}.nd2 = unique(ids(cid, 2));
				c{nc}.nn2 = length(c{nc}.nd2);
				c{nc}.id = sub2ind(sz, x1(cid), x2(cid), x3(cid));
			end
		end
	end
end

%% subfunction 4 edge clusters
function [c, nc] = cluster4Edge(in, nb)
	sz =size(in);
	[x1, x2, x3] = ind2sub(sz, find(in));
	nds = length(x1);
	fg = true(nds, 1);
	ids = [x1, x2, x3];
	nc = 0;
	c = [];
	for ig = 1 : nds
		if fg(ig)
			cid = [];
			[cid, fg] = searchClusterLoopSparseEdge(ids, ig, cid, fg, nb);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :);
				c{nc}.nd1 = unique(ids(cid, 1));
				c{nc}.nn1 = length(c{nc}.nd1);
				c{nc}.nd2 = unique(ids(cid, 2));
				c{nc}.nn2 = length(c{nc}.nd2);
				c{nc}.id = sub2ind(sz, ids(cid, 1), ids(cid, 2), ids(cid, 3));
			end
		end
	end
end

%% subfunction 5 time filter
function	nTime = timeFilter(in)
	% % % input size
	sz = size(in);	%expect ([n, n, T])
	tmp = zeros([sz, 2]);
	% % % shift array to 	
	tmp(:, :, 2 : end, 1) = in(:, :, 1 : end - 1);
	tmp(:, :, 1 : end - 1, 2) = in(:, :, 2 : end);
	nTime = sum(tmp, 4);
end

%% subfunction 6 space filter
function	nSpace = spaceFilter(in, nb)

	% % % iput size
	sz = size(in);	%expect [n n T]
	nSpace = zeros(sz);
	% % % check consistency
	if sz(1) ~= length(nb)
		error('inconsistent grid number in connection and neighbor definition!');
	end
	% % % get neighbors
	for iNode = 1 : sz(1)
		idx = nb{iNode};	%all neighbors
		nSpace(iNode, :, :) = nSpace(iNode, :, :) + sum(in(idx, :, :), 1);
		nSpace(:, iNode, :) = nSpace(:, iNode, :) + sum(in(:, idx, :), 2);
		% % % auto connection should be excluded from neighbors
		id0 = squeeze(logical(in(iNode, iNode, :)));
		nSpace(iNode, iNode, id0) = nSpace(iNode, iNode, id0) - 2;	
	end
end

%% subfunction 7 node finder
function [y, flag] = searchClusterLoopSparseNode(c3, x, y, flag)
	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			s = find(sum(bsxfun(@eq, c3, c3(x(k), :)), 2) >= 2);
			s(s==x(k)) = [];
			d = abs(bsxfun(@minus, c3(s, 3), c3(x(k), 3))) > 1;
			s(d) = [];
			x2 = [x2; s];
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparseNode(c3, x2, y, flag);
	end
end

%% subfunction 8 Edge finder
function [y, flag] = searchClusterLoopSparseEdge(c3, x, y, flag, nb)

	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			s1 = find(bsxfun(@eq, c3(:, 3), c3(x(k), 3)));
			if ~isempty(s1)
				for k1 = 1 : 2
					s2 = find(ismember(c3(s1, k1), c3(x(k), k1)));
					if ~isempty(s2)
						s3 = find(ismember(c3(s1(s2), 3 - k1), nb{c3(x(k), 3 - k1)}));
						if ~isempty(s3)
							x2 = [x2; s1(s2(s3))];
						end
					end
				end
			end
			clear s1 s2 s3 k1;
			
			s1 = find(sum(bsxfun(@eq, c3(:, 1 : 2), c3(x(k), 1 : 2)), 2) == 2);
			if ~isempty(s1)
				s2 = abs(bsxfun(@minus, c3(s1, 3), c3(x(k), 3))) == 1;
				if ~isempty(s2)
					x2 = [x2; s1(s2)];
				end
			end
			clear s1 s2;
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparseEdge(c3, x2, y, flag, nb);
	end
end

%% subfunction 9 time dimension summary
function [tmd, tmt] = cluster_sumTime(ed, nTs, t, id)
	% % % check inputs
	if nargin ~= 2 && nargin ~= 4 
		error('I need two inputs or four!');
	end

	% % % t-values flag
	if nargin < 3 || isempty(t)
		tFlag = false;
		tmt = [];
	else
		tFlag = true;
		tmt = zeros(nTs, 1);
	end
	
	% % % number of connections
	tmd = zeros(nTs, 1);
	for it = 1 : nTs
		idx = find(ed(:, 3) == it);
		tmd(it) = length(idx); 
		if tFlag
			tmt(it) = sum(t(id(idx)));
		end
	end
	
end %end of function

%% subfunction 10 spatial dimension summary
function [spd, spt] = cluster_sumSpace(ed, nvxs, t)
	% % % check inputs
	if nargin < 2 || nargin > 3
		error('I need two inputs or three!');
	end

	% % % t-values flag
	if nargin < 3 || isempty(t)
		tFlag = false;
		spt = [];
	else
		tFlag = true;
		spt = zeros(nvxs);
	end
	
	% % % number of connections
	spd = zeros(nvxs);
	for k =  1 : size(ed, 1)
		spd(ed(k, 1), ed(k, 2)) = spd(ed(k, 1), ed(k, 2)) + 1;
		if tFlag
			spt(ed(k, 1), ed(k, 2)) = spt(ed(k, 1), ed(k, 2)) + t(ed(k, 1), ed(k, 2), ed(k, 3));
		end
	end
end %end of function

%% subfunction 11 get the p evaluation
function p = cluster_pValue(H, sz)
	[y, I] = min(abs(H - sz));
	p = I ./ length(H);
	clear y I H sz;
end %end of function


%% subfunction 12 summarize all
function	c = sum_clusters(r, pp, t)
	% % % initiate
	cc = cell(1, 1);
	ct = 0;
	sz = zeros(1);
	
	% % % find
	for k = 1 : 2
		for ik = 1 : numel(r.c{k})
			if r.c{k}{ik}.sz > r.th
				ct = ct + 1;
				cc{ct} = r.c{k}{ik};
				cc{ct}.p0 = cluster_pValue(r.H, r.c{k}{ik}.sz);
				sz(ct) = cc{ct}.sz;
				if k > 1.5
					cc{ct}.sgn = -1;
				else
					cc{ct}.sgn = 1;
				end
				[cc{ct}.tmd, cc{ct}.tmt] = cluster_sumTime(cc{ct}.ed, pp.nTs, t, cc{ct}.id);
				[cc{ct}.spd, cc{ct}.spt] = cluster_sumSpace(cc{ct}.ed, pp.n, t);
			end
		end
	end
	
	% % %  sort the clusters
	nc = length(cc);
	if nc > 1
		[tmp, I] = sort(sz, 'descend');
		c = cc(I);
		clear tmp I cc;
	else
		c = cc;
	end
	
	% % % clean up
	clear t r pp sz;
end