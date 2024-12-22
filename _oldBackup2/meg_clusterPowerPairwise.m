function [c, th, H] = meg_clusterPowerPairwise(data, nb, th1, th2, nClustersPerTest, nTestPerCluster, nblimit, pth, fname, pFlag)

% % % rewritten 30/01/2019 for pairwise cluster
% % % updated 29/01/2019 move summary part out as stand-alone function
% % % updated 21/01/2019 to solve emepty clusters
% % % updated 20/12/2017 to new naming rules
% % % updated 07/02/2017 to allow bigger subject number (up to 90)
	

	%% data preparison
	% % % para defaults
	if nargin < 2
		error('At least four inputs needed: connection array, neigbor definition.');
	else
		% % % data [f, t, ch, sb]
		[nFs, nTs, nvxs, nSubs, pair] = size(data);
		if pair ~= 2
			fprintf('We need a pair of data!');
			c =[];th=[];H=[];
			return;
		end
	end
	if nargin < 3 || isempty(th1)
		th1 = 2.86;
	elseif th1 < 1
		th1 = mat_tThreshFromP(nSubs, th1);
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
	
	if nargin < 10 || isempty(pFlag)
		pFlag = 1;
	end
	
	%% initialize
	H0 = zeros(nClustersPerTest, nTestPerCluster);
	if nTestPerCluster > 2^nSubs
		error('The number of subjects is too small!');
	elseif nSubs <= 20
		s0 = uint32(randperm(2.^nSubs));
		s1 = ismember(dec2bin(s0(1 : nTestPerCluster) - 1, nSubs), '1');
		clear s0;
	else
		s1 = rand(nTestPerCluster, nSubs) > 0.5;
	end
	
	%% do computation
	t = mat_tValue4Arrays(data, 4, 0);
	[myClusters1, sz1, nc1] = theCluster4Power (t, nb, th1, th2, nblimit, pFlag);
% 	c = myClusters1{1};return; %%% debug only
	[myClusters2, sz2, nc2] = theCluster4Power (-t, nb, th1, th2, nblimit, pFlag);
	save(fname, 'myClusters1', 'sz1', 'nc1', 'myClusters2', 'sz2', 'nc2');
	
	% % % we took t without the control set	
	t = t(:, :, :, 1);
	
	%% get statistics
	parfor ts = 1 : nTestPerCluster
% 		for ts = 1 : nTestPerCluster
		tic;

		% % % prepare data		
		cz = bsxfun(@times, data, 2 * (permute(s1(ts, :), [1 3 4 2]) - 0.5));
		con = mat_tValue4Arrays(cz, 4, 0);

		% % % do cluster
		[tmp1, sz0, tmp3] = theCluster4Power (con, nb, th1, th2, nblimit, pFlag);
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
	[c, th, H] = meg_clusterSummaryPower(H0, myClusters1, myClusters2, nTs, nFs, nvxs, pth, t);	
	save(fname, 'c', 'th', 'H', '-append');
	
end %end of function



%% below are sub-functions needed
%% 1.main
function [myClusters, sz, nc] = theCluster4Power (conArray, neighborDefinition, conThresh, spfThresh, nblimit, pFlag)
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
	switch pFlag
		case 1 %network 1 only
			b1 = conArray(:, : , :, 1) > conThresh;
			b2 = conArray(:, : , :, 2) > conThresh;
			binArray = (b1 & ~b2);
		case 2 %network 1 only or network 1 is larger than network 2
			b1 = conArray(:, : , :, 1) > conThresh;
			b2 = conArray(:, : , :, 2) > conThresh;
			b3 = conArray(:, : , :, 1) - conArray(:, : , :, 2) > conThresh;
			binArray = (b1 & ~b2) | b3;
		case 3 %network 1 is larger than network 2
			binArray = conArray(:, : , :, 1) - conArray(:, : , :, 2) > conThresh;
		otherwise
			binArray = conArray(:, : , :, 1) > conThresh;
	end
			
	clear b1 b2 b3;
	
	%% spatial filtering
	binArray = spatialFilter4PowerCluster(binArray, neighborDefinition, spfThresh);
	
	%% find cluster
	if nblimit > 0
		[myClusters, nc] = findCluster4PowerNeighborRestrainedPowerCluster (binArray, neighborDefinition);
	elseif nblimit < 0
		[myClusters, nc] = findCluster4PowerNeighborStronglyRestrainedPowerCluster (binArray, neighborDefinition, -nblimit);
	else
		[myClusters, nc] = findCluster4PowerNeighborFree (binArray);
	end
	
	%% sort cluster
	if nc
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
		myClusters = [];
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
	tmpArray(2 : end, :, :, 1) = binArray(1 : end - 1, :, :);
	tmpArray(1 : end - 1, :, :, 2) = binArray(2 : end, :, :);
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



%% 15. nb restrained cluster
function [c, nc] = findCluster4PowerNeighborStronglyRestrainedPowerCluster (binArray, nb, th)

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
			[cid, flag] = searchClusterLoopSparser(ids, ig, cid, flag, nb, th);
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

%% 16. nb restrained search
function [y, flag] = searchClusterLoopSparser(ids, x, y, flag, nb, th)
	NBSHIFT = [0 0 1; 0 1 0; 0 0 -1; 0 -1 0];
	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			% % % spatial neighbor, same as before
			s1 = ismember(ids(:, 2:3), ids(x(k), 2:3), 'rows');
			s2 = ismember(ids(:, 1), nb{ids(x(k), 1)});
			x2 = [x2; find(s1 & s2)];
			clear s1 s2;
			% % % time-frequency neighbor, more restrict
			for ic = 1 : 4
				s1 = ids(x(k), :) + NBSHIFT(ic, :);
				[s5, s6] = ismember(s1, ids, 'rows');
				if s5
					s2 = repmat(nb{ids(x(k), 1)}, [1 3]);
					s2(:, 2) = s1(2);
					s2(:, 3) = s1(3);
					s4 = ismember(s2, ids, 'rows');
					s7 = numel(s4(s4)) / numel(s4);
					if s7 > th
						x2 = [x2; s6];
					end
				end
				clear s1 s2 s3 s4 s5 s6 s7;
			end
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparser(ids, x2, y, flag, nb, th);
	end
end

