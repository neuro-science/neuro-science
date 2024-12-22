function [c, th, H] = meg_clusterPowerSingle(data, trls, nb, th1, th2, nClustersPerTest, nTestPerCluster, nblimit, pth, fname)

% % % updated 20/12/2017 to new naming rules
% % % updated 07/02/2017 to allow bigger subject number (up to 90)
	

	%% data preparison
	% % % para defaults
	if nargin < 2
		error('At least four inputs needed: connection array, neigbor definition.');
	else
		% % % data [f, t, ch, sb]
		[nFs, nTs, nvxs, nTrls] = size(data);
	end
	if nargin < 3 || isempty(th1)
		th1 = 2.86;
	elseif th1 < 1
		th1 = mat_tThreshFromP(nTrls, th1);
	end
	if nargin < 4 || isempty(th2)
		th2 = 0.5;
	end
	if nargin < 5 || isempty(nClustersPerTest)
		nClustersPerTest = 3;
	end
	if nargin < 6 || isempty(nTestPerCluster)
		nTestPerCluster = 1000;
	end
	if nargin < 7 || isempty(nblimit)
		nblimit = 1;
	end
	
	if nargin < 8 || isempty(pth)
		pth = 0.05;
	end
	
	if nargin < 9 || isempty(fname)
		fname = 'unknown_cluster4pwr.mat';
	end
	
	%% initialize
	H0 = zeros(nClustersPerTest, nTestPerCluster);
	s1 = rand(nTestPerCluster, nTrls) > 0.5;
	
	%% do computation
	t = mat_tValue4ArraysIndependent(data, trls, 4);
	t = reshape(t, [nFs, nTs, nvxs]);
	[myClusters1, sz1, nc1] = theCluster4Power (t, nb, th1, th2, nblimit);
	[myClusters2, sz2, nc2] = theCluster4Power (-t, nb, th1, th2, nblimit);
	save(fname, 'myClusters1', 'sz1', 'nc1', 'myClusters2', 'sz2', 'nc2');
	
	%% get statistics
	for ts = 1 : nTestPerCluster
% 		for ts = 1 : nTestPerCluster
		tic;

		% % % prepare data		
		con = mat_tValue4ArraysIndependent(data, s1(ts, :), 4);
		con = reshape(con, [nFs, nTs, nvxs]);

		% % % do cluster
		[tmp1, sz0, tmp3] = theCluster4Power (con, nb, th1, th2, nblimit);
		tmp1 = [];
		tmp3 = [];

		% % % fill the results		
		szn = zeros(nClustersPerTest, 1);
		ssz = length(sz0);
		if ssz >= nClustersPerTest
			szn(:, 1) = sz0(1 : nClustersPerTest);
		elseif ssz > 0
			szn(1 : ssz, 1) = sz0;
		end
		H0(:, ts) = szn;

		% % % echo message		
		fprintf('Test #%03d done after %7.2f seconds.\n', ts, toc);
	end
	save(fname, 'H0', '-append');
	
	%% summary
	[c, th, H] = clusterSummaryPower(H0, myClusters1, myClusters2, nTs, nFs, nvxs, pth, t);	
	save(fname, 'c', 'th', 'H', '-append');
	
end %end of function



%% below are sub-functions needed
%% 1.main
function [myClusters, sz, nc] = theCluster4Power (conArray, neighborDefinition, conThresh, spfThresh, nblimit)
	%cA[f, t, ch]
	%% preparision
	if nargin < 2
		error('At least four inputs needed: connection array, neigbor definition.');
	elseif nargin < 3
		conThresh = 2.86;
	elseif nargin < 4
		spfThresh = 0.5;
	elseif nargin < 5
		nblimit = 1;
	end
	
	%% thresholding to binary
	binArray = conArray > conThresh;
	
	%% spatial filtering
	binArray = spatialFilter4PowerCluster(binArray, neighborDefinition, spfThresh);
	
	%% find cluster
	if nblimit
		[myClusters, nc] = findCluster4PowerNeighborRestrainedPowerCluster (binArray, neighborDefinition);
	else
		[myClusters, nc] = findCluster4PowerNeighborFree (binArray);
	end
	
	%% sort cluster
	if nc > 0
		tmp = zeros(nc, 1);
		nn = zeros(nc, 1);
		for k = 1 : nc
			tmp(k) = sum(conArray(myClusters{k}.id));
			myClusters{k}.sz = tmp(k);
			nn(k) = myClusters{k}.nn;
		end
		[sz, I] = sort(tmp, 'descend');
		myClusters = myClusters(I);
		clear I tmp;
	else
		sz = [];
	end
end

%% 2. spatial filter
function	[outBinArray, nNeighborsInTimeFreq, nNeighborsTopoSpatial] = ...
	spatialFilter4PowerCluster(binArray, neighborDefinition, spfThresh)
% % % modified 07/02/2017 to fit 1t1f condition - to be done
% % % modified 4/11/2013 both end and 1st

	% % % prepare para
	sz = size(binArray);	%expect(f, t, node)
	nNodes = length(neighborDefinition);
	
	% % % numbers of possible spatial neighbors	
	nNeighborsPossiblePerNode = zeros(nNodes, 1);
	for k = 1 : nNodes
		nNeighborsPossiblePerNode(k) = length(neighborDefinition{k}) - 1;	%spatial neighbor per node
	end
	
	% % % numbers of possible time-frequency neighbors	
	nNeighborsPossiblePerTimeFreq = zeros(sz(1 : 2)) + 4;
	nNeighborsPossiblePerTimeFreq(1, :) = nNeighborsPossiblePerTimeFreq(1, :) - 1;
	nNeighborsPossiblePerTimeFreq(:, 1) = nNeighborsPossiblePerTimeFreq(:, 1) - 1;
	nNeighborsPossiblePerTimeFreq(end, :) = nNeighborsPossiblePerTimeFreq(end, :) - 1;
	nNeighborsPossiblePerTimeFreq(:, end) = nNeighborsPossiblePerTimeFreq(:, end) - 1;

	% % % real number of Neighbors	
	nNeighborsInTimeFreq = tfFilter(binArray);
	nNeighborsTopoSpatial = topoFilter(binArray, neighborDefinition);
	
	% % % sum together and apply filter	
	nNeighborsPossibleAll = bsxfun(@plus, permute(nNeighborsPossiblePerNode, [2 3 1]), nNeighborsPossiblePerTimeFreq);
	outBinArray = binArray & ((nNeighborsInTimeFreq + nNeighborsTopoSpatial) > spfThresh * nNeighborsPossibleAll);
end

%% 3. tf filter
function nNeighborsInTimeFreq = tfFilter(binArray)
	sz = size(binArray);	%expect ([f, t, nodes)
	tmpArray = zeros([sz, 4]);
	% % % shift array to 	
	tmpArray(2 : end, :, :, 1) = binArray(1 : end - 1, : ,:);
	tmpArray(1 : end - 1, :, :, 2) = binArray(2 : end, : ,:);
	tmpArray(:, 2 : end, :, 3) = binArray(:, 1 : end - 1, :);
	tmpArray(:, 1 : end - 1, :, 4) = binArray(:, 2 : end, :);
	nNeighborsInTimeFreq = sum(tmpArray, 4);
end

%% 4. topo filter
function	nNeighborsTopoSpatial = topoFilter(binArray, neighborDefinition)

	sz = size(binArray);	%expect(f, t, node)
	nNeighborsTopoSpatial = zeros(sz);
	if sz(3) ~= length(neighborDefinition)
		error('inconsistent grid number in connection and neighbor definition!');
	end
	
	for iNode = 1 : sz(3)
		idx = neighborDefinition{iNode};	%all neighbors
		idx(idx == iNode) = [];	%exclude itself
		nNeighborsTopoSpatial(:, :, iNode) = nNeighborsTopoSpatial(:, :, iNode) + sum(binArray(:, :, idx), 3);
	end
% 	nNeighborsTopoSpatial(nNeighborsTopoSpatial < 0) = 0;
end

%% 5. nb restrained cluster
function [c, nc] = findCluster4PowerNeighborRestrainedPowerCluster (binArray, nb)

	sz =size(binArray);
	[x1, x2, x3] = ind2sub(sz, find(binArray));
	nds = length(x1);
	flag = true(nds, 1);
	ids = [x3, x1, x2];
	nc = 0;
	c = [];
	for ig = 1 : nds
		if flag(ig)
			cid = [];
			[cid, flag] = searchClusterLoopSparse(ids, ig, cid, flag, nb);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :);
				c{nc}.nd = unique(ids(cid, 1));
				c{nc}.nn = length(c{nc}.nd);
				c{nc}.id = sub2ind(sz, x1(cid), x2(cid), x3(cid));
			end
		end
	end
end

%% 6. nb restrained search
function [y, flag] = searchClusterLoopSparse(ids, x, y, flag, nb)
	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			d = sum(abs(bsxfun(@minus, ids(:, 2:3), ids(x(k), 2:3))), 2);
			s1 = ids(:, 1) == ids(x(k), 1);
			s2 = d == 1;
			s3 = d == 0;
			s4 = ismember(ids(:, 1), nb{ids(x(k), 1)});
			s = find((s1 & s2) | ((~s1) & s3 & s4));
			x2 = [x2; s];
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparse(ids, x2, y, flag, nb);
	end
end

%% 7. nb free cluster
function [c, nc] = findCluster4PowerNeighborFree (binArray)

	sz =size(binArray);
	[x1, x2, x3] = ind2sub(sz, find(binArray));
	nds = length(x1);
	flag = true(nds, 1);
	ids = [x3, x1, x2];
	nc = 0;
	c = [];
	for ig = 1 : nds
		if flag(ig)
			cid = [];
			[cid, flag] = searchClusterLoopSparseNoNB(ids, ig, cid, flag);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :);
				c{nc}.nd = unique(ids(cid, 1));
				c{nc}.nn = length(c{nc}.nd);
				c{nc}.id = sub2ind(sz, x1(cid), x2(cid), x3(cid));
			end
		end
	end
end

%% 8. nb free search
function [y, flag] = searchClusterLoopSparseNoNB(ids, x, y, flag)
	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			s = find(sum(bsxfun(@eq, ids, ids(x(k), :)), 2) >= 2);
			s(s==x(k)) = [];
			d = sum(abs(bsxfun(@minus, ids(s, 2:3), ids(x(k), 2:3))), 2) > 1;
			s(d) = [];
			x2 = [x2; s];
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparseNoNB(ids, x2, y, flag);
	end
end

%% 9. summary
function [cc, th, H] = clusterSummaryPower(H0, c1, c2, nTs, nFs, nvxs, pth, t)
% % % 15/07/2015	updated by wp, t
% % % 26/02/2015	written by wp

	%% prepare
	% % % paras
	if nargin < 8 || isempty(t)
		tflag = false;
	elseif any([size(t, 1) - nFs, size(t, 2) - nTs, size(t, 3) - nvxs])
		str1 = sprintf('%d ', size(t));
		str2 = sprintf(' -vs.- ');
		str3 = sprintf('%d ', [nFs, nTs, nvxs]);
		str4 = sprintf(': wrong size of t!\n');
		str = [str1, str2, str3, str4];
		error(str);
	else
		tflag = true;
	end
	
	% % % paras
	if nargin < 7 || isempty(pth)
		pth = 0.05;
	end
	
	if nargin < 6 || isempty(nvxs)
		nvxs = 324;
	end
	
	if nargin < 5 || isempty(nFs)
		nFs = 21;
	end
	
	if nargin < 4 || isempty(nTs)
		nTs = 17;
	end
	
	if nargin < 3
		c2 = [];
	end
	
	if nargin < 2
		error('At leaset 2 inputs are needed!');
	end
	
	% % % data
	c0{1} = c1;
	c0{2} = c2;
	
	%% post-process
	% % % get the thresh
	[th, H] = mat_distr2thresh(H0, pth);
	
	% % % select the clusters
	cc = cell(1, 1);
	ct = 0;
	for k = 1 : 2
		for ik = 1 : numel(c0{k})
			if c0{k}{ik}.sz > th
				ct = ct + 1;
				cc{ct} = c0{k}{ik};
				cc{ct}.p0 = cluster_pValue(H, c0{k}{ik}.sz);
				if k > 1.5
					cc{ct}.sgn = -1;
				else
					cc{ct}.sgn = 1;
				end
				cc{ct}.tfd = cluster_TF2(cc{ct}.ed, nFs, nTs);
				cc{ct}.spd1 = cluster_SP1(cc{ct}.ed, nvxs);
				if tflag
					cc{ct}.spt1 = cluster_SPT1(cc{ct}.ed, nvxs, t);
					cc{ct}.tft = cluster_TFT2(cc{ct}.ed, nFs, nTs, t);
				end
			end
		end
	end
	
	% % %  sort the clusters
	nc = length(cc);
	if nc > 1
		sz1 = zeros(nc, 1);
		for k = 1 : nc
			sz1(k) = cc{k}.sz;
		end
		[tmp, I] = sort(sz1, 'descend'); %modeified 09/05/2017
		cc = cc(I);
		clear tmp I sz1;
	end
	
	fprintf('=========Data sorted on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));
end % end of function

%% 10. tf2d with n
function tfd = cluster_TF2(ed, nFs, nTs)
	tfd = zeros(nFs, nTs);
	for it = 1 : nTs
		for iq = 1 : nFs
			tfd(iq, it) = length(find(ed(:, 2) == iq & ed(:, 3) == it)); 
		end
	end
end %end of function

%% 11. TF2d with t
function tft = cluster_TFT2(ed, nFs, nTs, t)
	tft = zeros(nFs, nTs);
	for k = 1 : size(ed, 1)
		tft(ed(k, 2), ed(k, 3)) = tft(ed(k, 2), ed(k, 3)) + t(ed(k, 2), ed(k, 3), ed(k, 1)); 
	end
end %end of function

%% 12. sp 1d with n
function spd1 = cluster_SP1(ed, nvxs)
	tmp = ed(:, 1);
	spd1 = zeros(nvxs, 1);
	for k =  1 : size(tmp, 1)
		spd1(tmp(k)) = spd1(tmp(k)) + 1;
	end
	clear tmp;
end %end of function

%% 13. sp 1d with t
function spt1 = cluster_SPT1(ed, nvxs, t)
	spt1 = zeros(nvxs, 1);
	for k =  1 : size(ed, 1)
		spt1(ed(k, 1)) = spt1(ed(k, 1)) + t(ed(k, 2), ed(k, 3), ed(k, 1));
	end
	clear tmp;
end %end of function

%% 14. cluster p value
function p = cluster_pValue(H, sz)
	[y, I] = min(abs(H - sz));
	p = I ./ length(H);
	clear y I H sz;
end %end of function


