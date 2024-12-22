function c = meg_clusterPowerANOVA(data, nb, th1, th2, nClustersPerTest, nTestPerCluster, nblimit, pth, fname)
% % % rewritten on 27/02/2019, for Fthresh
% % % rewritten on 06/02/2019, for better anova
% % % rewritten on 04/02/2019, for summary
% % % rewritten on 30/01/2019, for newer naming rules
% % % written on 08/05/2017, based on the cmp30_Permute4H0_Cluster4Power
	

	%% 1. data preparison
	%% data preparison
	cstr = {'F1Main', 'F2Main', 'F1*F2I'};
	% % % para defaults
	if nargin < 2
		error('At least two inputs needed: connection array, neigbor definition.');
	else
		% % % data [f, t, ch, sb, lf1, lf2], here lf1 and lf2 are the number of levels for factor 1 and 2
		[nFs, nTs, nvxs, nSubs, lf1, lf2] = size(data);	
	end
	
	if nargin < 3 || isempty(th1)
		th1F = repmat(2.86, [3 1]);
	elseif numel(th1) == 1
		th1 = repmat(th1, [3, 1]);
	elseif numel(th1) ~= 3
		fprintf('wrong size of thresh 1, exiting...\n');
	end
		
	if max(th1) > 1
		th1F = th1;
	else
		nn{1} = lf1 - 1;
		nn{2} = lf2 - 1;
		nn{3} = nn{1} * nn{2};
		for ic = 3 : -1 : 1
			dn = (nSubs - 1) * nn{ic};
			th1F(ic) = mat_FThreshFromP(nn{ic}, dn, th1(ic));
			fprintf('I compute F threshold %s from P-values: th = %7.3f, as p = %7.5f, N = %d.\n', ...
				cstr{ic}, th1F(ic), th1(ic), nSubs);
			clear dn;
		end
		clear nn;
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
		
	
	%% 2. do computation
	% % % the real data
	F = rmANOVA4ThisData(data);	% % % size: (nFs, nTs, nvxs, 3);
	for ic = 3 : -1 : 1
		[r.cl{ic}, r.sz{ic}, r.nc{ic}] = theCluster4Power (F(:, :, :, ic), nb, th1F(ic), th2, nblimit);
	end
	save(fname, 'r', 'F');
	% % % statistics of null distribution
	H0 = zeros(nClustersPerTest, nTestPerCluster, 3);
	if nTestPerCluster > 2^nSubs
		error('The number of subjects is too small!');
	elseif nSubs <= 12
		rng('shuffle');
		s0 = uint32(randperm(2.^nSubs));
		s1(:, :, 2) = ismember(dec2bin(s0(1 : nTestPerCluster) - 1, nSubs), '1');
		rng('shuffle');
		s0 = uint32(randperm(2.^nSubs));
		s1(:, :, 1) = ismember(dec2bin(s0(1 : nTestPerCluster) - 1, nSubs), '1');
		clear s0;
	else
		rng('shuffle');
		s1 = rand(nTestPerCluster, nSubs, 2) > 0.5;
	end
	s1 = permute(s1, [2 3 1]); %[test, sub, 2] - [sub, 2, test]

	parfor ts = 1 : nTestPerCluster
%	for ts = 1 : nTestPerCluster
		tic;
		% % % prepare data		
		tF = rmANOVA4ThisData(data, s1(:, :, ts));	% % % size: (nFs, nTs, nvxs, 3);

		for ic = 1 : 3
			% % % do cluster
			[tmp1, sz0, tmp3] = theCluster4Power (tF(:, :, :, ic), nb, th1F(ic), th2, nblimit);
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
			H0(:, ts, ic) = szn;

			% % % echo message		
			fprintf('Test #%03d done after %7.2f seconds.\n', ts, toc);
		end
	end
	save(fname, 'H0', '-append');
	
	%% 3. summarize the data
	% % % define some parameters
	cstr = {'F1main', 'F2main', 'F1*F2'};
	% % % initialize
	th = zeros(3, 1);
	cc = cell(2, 1);
	ct = 0;
	% % % gether the data in a loop
	for ic = 1 : 3
		% % % get the thresh
		[th(ic), H] = mat_distr2thresh(H0(:, :, ic), pth);
		% % % check the clusters one by one
		for ik = 1 : numel(r.cl{ic})
			if r.cl{ic}{ik}.sz > th(ic)
				ct = ct + 1;
				cc{ct} = r.cl{ic}{ik};
				cc{ct}.prefix = cstr{ic};
				cc{ct}.p0 = cluster_pValue(H, r.cl{ic}{ik}.sz);
				if nTs * nFs > 1
					cc{ct}.tfd = cluster_TF2(cc{ct}.ed, nFs, nTs);
					cc{ct}.tff = cluster_TFT2(cc{ct}.ed, nFs, nTs, F(:, :, :, ic));
				end
				cc{ct}.spd1 = cluster_SP1(cc{ct}.ed, nvxs);
				cc{ct}.spf1 = cluster_SPT1(cc{ct}.ed, nvxs, F(:, :, :, ic));
			end
		end
	end
	
	% % %  sort the clusters
	if ct
		sz1 = zeros(ct, 1);
		for k = 1 : ct
			sz1(k) = cc{k}.sz;
		end
		[tmp, I] = sort(sz1, 'descend');
		c = cc(I);
		clear tmp I sz1 cc;
	else
		c = [];
	end
	
	% % % save and notify
	save(fname, 'c', 'th', '-append');
	fprintf('=========Data sorted on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));	
end %end of function

function F = rmANOVA4ThisData(data, rFlag)
	
	d = permute(double(data), [5 6 4 1 2 3]);	%[f, t, ch, sb, lf1, lf2] -> [lf1, lf2, sb, f, t, ch]
	if nargin < 2 || isempty(rFlag) 
		F = mat_rmANOVA2Array(d);	%	[f, t, ch, 3]
	else	% permute the data for control distribution
		d(:, :, rFlag(:, 1), :, :, :) = d([2 1], :, rFlag(:, 1), :, :, :);
		d(:, :, rFlag(:, 2), :, :, :) = d(:, [2 1], rFlag(:, 2), :, :, :);
		F = mat_rmANOVA2Array(d);	%	[f, t, ch, 3]
	end
	% % % 	%% 1. check data and reformat
	% % % 	% % % data [f, t, ch, sb, lf1, lf2], here lf1 and lf2 are the number of levels for factor 1 and 2
	% % % 	[nFs, nTs, nvxs, nSubs, lf1, lf2] = size(data);	
	% % % 	% % % reshape the data to meet the requirement
	% % % 	data = reshape(data, [nFs * nTs * nvxs, nSubs* lf1 * lf2])'; %[num of voxels, num of tests]
	% % % 	% % % parameters for the case-definition
	% % % 	sb = repmat(1 : nSubs, [1, lf1 * lf2])'; %subjects
	% % % 	if rFlag % permute the data for control distribution
	% % % 		seq = zeros(nSubs, lf1*lf2);
	% % % 		for k = 1 : nSubs
	% % % 			seq(k, :) = randperm(lf1*lf2);
	% % % 		end
	% % % 	else % the real data
	% % % 		seq = repmat(1 : lf1*lf2, [nSubs, 1]);
	% % % 	end
	% % % 	f1 = mod(seq - 1, lf1) + 1; f1 = f1(:); %factor 1
	% % % 	f2 = ceil(seq / lf1); f2 = f2(:);	%factor 2
	% % % 	
	% % % 	%% 2. compute data
	% % % 	% % % initialize
	% % % 	F = zeros(nFs * nTs * nvxs, 3);
	% % % 	p = zeros(nFs * nTs * nvxs, 3);
	% % % 	% % % compute: a lot of loops, better use parallel computing
	% % % 	if multiCoreFlag
	% % % 		parfor iv = 1 : nFs * nTs * nvxs
	% % % 			stats = mat_rmANOVA2(data(:, iv), sb, f1, f2, {'F1', 'F2'});
	% % % 			F(iv, :) = [stats{2, 5}, stats{3, 5}, stats{4, 5}];
	% % % 			p(iv, :) = [stats{2, 6}, stats{3, 6}, stats{4, 6}];
	% % % 		end
	% % % 	else
	% % % 		for iv = 1 : nFs * nTs * nvxs
	% % % 			stats = mat_rmANOVA2(data(:, iv), sb, f1, f2, {'F1', 'F2'});
	% % % 			F(iv, :) = [stats{2, 5}, stats{3, 5}, stats{4, 5}];
	% % % 			p(iv, :) = [stats{2, 6}, stats{3, 6}, stats{4, 6}];
	% % % 		end
	% % % 	end
	% % % 	
	% % % 	%% 3. output: reshape the results
	% % % 	F = reshape(F, nFs, nTs, nvxs, 3);
	% % % 	p = reshape(p, nFs, nTs, nvxs, 3);
end


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
