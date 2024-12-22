%% function #1
function c = meg_clusterCoupling(data, nb, ...
	th1, th2, nClustersPerTest, nTestPerCluster, method, ...
	oid, npools, fileName, nTs, nFs, nVs, th3)
	% % % updated 26/02/2019 by wp : improve randomization
	% % % updated 02/07/2018 by wp : small message bug corrected
	% % % updated 14/11/2017 by wp : merge several sub functions
	% % % updated 09/03/2017 by wp : add parallel computation
	% % % updated 26/11/2015 by wp : add output 

	%% 01. check input and set defaults
	% % % para defaults
	if nargin < 2 || isempty(data) || isempty(nb)
		error('At least two inputs needed: connection array, neigbor definition.');
	else
		% % % check data size [f, t, ch, ch, sb]
		[nfs, nts, nvs1, nvs2, nSubs] = size(data);
		if nvs1 == nvs2
			nvs = nvs1;
			clear nvs1 nvs2;
		else
			error('data size mismatch!');
		end
	end
	if nargin < 3 || isempty(th1)
		clp.th1 = 0.05;
	else
		clp.th1 = th1;
	end
	if nargin < 4 || isempty(th2)
		clp.th2 = 0.5;
	else
		clp.th2 = th2;
	end
	if nargin < 5 || isempty(nClustersPerTest)
		clp.nClustersPerTest = 3;
	else
		clp.nClustersPerTest = nClustersPerTest;
	end
	if nargin < 6 || isempty(nTestPerCluster)
		clp.nTestPerCluster = 1000;
	else
		clp.nTestPerCluster = nTestPerCluster;
	end
	if nargin < 7 || isempty(method)
		clp.method = 'node';
	else
		clp.method = method;
	end
	if nargin < 8 || isempty(oid)
		oid = 1;
	end
	if nargin < 9 || isempty(npools)
		clp.npools = 1;
	else
		clp.npools = npools;
	end
	if nargin < 10 || isempty(fileName)
		clp.fileName = sprintf('tmpClusterCoupling_CTH%03dSTH%02d_%s_%04d-%02d-%02d_%02d_%02d_%02d.mat', clp.th1*1000, clp.th2*100, clp.method, round(clock));
	else
		clp.fileName = fileName;
	end
	if nargin < 11 || isempty(nTs)
		clp.nTs = nts;
	else
		clp.nTs = nTs;
	end
	if nargin < 12 || isempty(nFs)
		clp.nFs = nfs;
	end
	if nargin < 13 || isempty(nVs)
		clp.nVs = nvs;
	else
		clp.nVs = nVs;
	end
	if nargin < 14 || isempty(th3)
		clp.th3 = 0.05;
	else
		clp.th3 = th3;
	end
	
	%% 02. data preparison
	% % % connection threshold
	if clp.th1 > 1
		fprintf('\n\n!!!Warning, the connection threshold is larger than 1, I assume legacy values@%7.3f!\n\n', clp.th1);
		clp.th1t = clp.th1;
	else
		clp.th1t = mat_tThreshFromP(nSubs, clp.th1);
		fprintf('I compute threshold from P-values: th = %7.3f, as p = %7.5f, N = %d.\n', clp.th1t, clp.th1, nSubs);
	end
	% % % initialize H0
	H0 = zeros(clp.nClustersPerTest, clp.nTestPerCluster);
	% % % generate randomization
	tic;
	rng('shuffle');
	if nSubs < 12.5
		s0 = randperm(2.^nSubs - 2);
		s1 = ismember(dec2bin(s0(1 : clp.nTestPerCluster), nSubs), '1');
		clear s0;
	else
		s1 = rand(clp.nTestPerCluster, nSubs) > 0.5;
	end
	fprintf('Shuffle cost %6.1f seconds!\n', toc);
	% % % save the data for now
	save(clp.fileName, 'clp');
	
	%% 03. work on real data
	% % % connection real t-value
	con = mat_tValue4Arrays(data, 5, 0); %[f, t, ch, ch]
	% % % work on positive change
	fprintf(oid, '\n=========%s: Data computed in on %04d-%02d-%02d %02d:%02d:%02d=========\n', clp.method, round(clock));
	[myClusters1, sz1, nc1] = cluster4Coupling (con, nb, clp.th1t, clp.th2, clp.method);
	% % % work on negative change
	[myClusters2, sz2, nc2] = cluster4Coupling (-con, nb, clp.th1t, clp.th2, clp.method);
	fprintf(oid, '\n=========%s: Data computed out on %04d-%02d-%02d %02d:%02d:%02d=========\n', clp.method, round(clock));
	% % % save the data for now
	save(clp.fileName, 'myClusters1', 'sz1', 'nc1', 'myClusters2', 'sz2', 'nc2', '-append');
	% % % return if no clusters found
	if isempty(myClusters1) && isempty(myClusters2)
		fprintf('\n\n\nNo clusters found!!!\n\n\n');
		c = [];
		return;
	else
		fprintf('\nThe size: %10.2f!\n', max(abs([sz1(:); sz2(:)])));
	end
	
	%% 04. work on distribution of null data in loop
	if clp.npools > 1 % multiple CPU available
		% % % clean up cpus
		poolobj = gcp('nocreate');
		delete(poolobj);	
		% % % start the cpus
		thePool = parpool(clp.npools); 
		% % % run in parallel
		parfor ts = 1 : clp.nTestPerCluster
			tic;
			% % % prepare data		
			cz = bsxfun(@times, data, 2 * (permute(s1(ts, :), [1 3 4 5 2]) - 0.5));
			con_r = mat_tValue4Arrays(cz, 5, 0);
			cz = [];

			% % % do cluster
			[tmp1, sz0, tmp2] = cluster4Coupling (con_r, nb, clp.th1t, clp.th2, clp.method);
			tmp1 = [];
			tmp2 = [];

			% % % check the results and give feedback		
			szn = zeros(clp.nClustersPerTest, 1);
			ssz = length(sz0);
			if ssz >= clp.nClustersPerTest
				szn(:, 1) = sz0(1 : clp.nClustersPerTest);
				fprintf(oid, 'Hx = %10.1f.\nTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
			elseif ssz > 0
				szn(1 : ssz, 1) = sz0;
				fprintf(oid, 'Hx = %10.1f.\nTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
			else
				fprintf(oid, 'Hx = NaN.\nTest #%03d done after %7.2f seconds.\n', ts, toc);
			end
			
			% % % fill H0
			H0(:, ts) = szn;

		end %end of all repeats
		% % % clean up cpus
		delete(thePool);
	else
		for ts = 1 : clp.nTestPerCluster
			tic;
			% % % prepare data		
			cz = bsxfun(@times, data, 2 * (permute(s1(ts, :), [1 3 4 5 2]) - 0.5));
			con_r = mat_tValue4Arrays(cz, 5, 0);
			cz = [];

			% % % do cluster
			[tmp1, sz0, tmp2] = cluster4Coupling (con_r, nb, clp.th1t, clp.th2, clp.method);
			tmp1 = [];
			tmp2 = [];

			% % % check the results and give feedback		
			szn = zeros(clp.nClustersPerTest, 1);
			ssz = length(sz0);
			if ssz >= clp.nClustersPerTest
				szn(:, 1) = sz0(1 : clp.nClustersPerTest);
				fprintf(oid, 'Hx = %10.1f.\nTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
			elseif ssz > 0
				szn(1 : ssz, 1) = sz0;
				fprintf(oid, 'Hx = %10.1f.\nTest #%03d done after %7.2f seconds.\n', sz0(1), ts, toc);
			else
				fprintf(oid, 'Hx = NaN.\nTest #%03d done after %7.2f seconds.\n', ts, toc);
			end
			
			% % % fill H0
			H0(:, ts) = szn;

		end
	end
	clear ts szn ssz tmp1 tmp2 s1 con_r;
	
	% % % save the data for now
	save(clp.fileName, 'H0', '-append');
	
	%% 05. summarize the data
	
	% % % pool the data
	c0{1} = myClusters1;
	c0{2} = myClusters2;
	
	% % % get the thresh
	[th, H] = distr2thresh4H0(H0, clp.th3);
	
	% % % merge the clusters
	c = cell(1, 1); %assume not empty
	ct = 0;
	for k = 1 : 2
		for ik = 1 : numel(c0{k})
			if c0{k}{ik}.sz > th
				ct = ct + 1;
				c{ct} = c0{k}{ik};
				c{ct}.p0 = cluster_pValue(H, c0{k}{ik}.sz);
				if k > 1.5
					c{ct}.sgn = -1;
				else
					c{ct}.sgn = 1;
				end
				c{ct}.tfd = cluster_TF2(c{ct}.ed, clp.nFs, clp.nTs);
				[c{ct}.spd1, c{ct}.spd2] = cluster_SP2(c{ct}.ed, clp.nVs);
				[c{ct}.spt2, c{ct}.spt1, c{ct}.tft, c{ct}.te] = ed2tClusterCoupling(c{ct}.ed, con, clp.nVs, clp.nFs, clp.nTs);
			end
		end
	end
	
	% % %  sort the clusters
	nc = length(c);
	if nc > 1
		sz1 = zeros(nc, 1);
		for k = 1 : nc
			sz1(k) = c{k}.sz;
		end
		[tmp, I] = sort(sz1, 'descend');
		c = c(I);
		clear tmp I sz1;
	elseif isempty(c{1})
		c = [];
	end
	fprintf('=========Data sorted on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));
	
	% % % save the data
	save(clp.fileName, 'c', 'H', 'th', '-append');
	
end %end of function

%% function #2
function [myClusters, sz, nc] = cluster4Coupling (conArray, neighborDefinition, conThresh, spfThresh, method)

% cluster4Coupling will find clusters based on the 
% connection array, which defines cnnection strength between locations
% in various time and frequency. The algirthm came from the 2011 Neuron
% paper by Hipp, J.F. Engel, A.K. and Siegel, M.

% Use as
%   myClusters = cluster4Coupling (conArray, neighborDefinition, conThresh, spfThresh, method)
% 
% Input:
%		conArray					- a numFreq x numTime x numGridPair array, defining the
%									connection strength between two locations at certain 
%									time and frequency
%		neighborDefinition	- a numGrid x 1 cell array, k-th element contains neigboring 
%									grids of grid k.
%		conThresh				- threshhold of connection strength, value in conArray 
%									above it will indicate a connection or else not. 
%		spfThresh				- spatial filter will be applied to remove spurious connection.
%		method               - method of clustering, can be 'node' or 'edge' 
% 
% 
% Output:
%		myClusters - a cell array, each element contains a cluster
% 
% 
% 
% Written by Peng Wang, Institute of Neurophysiology, UKE.
%
%   Modification History:
%   00/08/2013: First Formal Draft
%   00/10/2013: improved Draft to allow more strategies
%   00/05/2014: improved Draft for more compact sub function
%   17/09/2014: updated by wp for better compatibility
%   14/11/2017: updated by wp for better compatibility


	%% 11. data preparison 
	% % % thresholding to binary
	binArray = conArray > conThresh;
% 	fprintf('1st floor!\n');	%###
	% % % spatial filtering
	binArray = 	spatialFilterCluster4Coupling (binArray, neighborDefinition, spfThresh);
% 	fprintf('2nd floor!\n'); %###
	
	%% 12. find cluster
	switch lower(method(1))
		case 'n'
			[myClusters, nc] = findCluster4CouplingNode(binArray);
		case 'e'
			[myClusters, nc] = findCluster4CouplingEdge(binArray, neighborDefinition);
		otherwise
			error('Current supported methods are <node> and <edge> only!');
	end
	
	%% 13. sort cluster
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
		J = find(nn > 10);
		clear I tmp;
	else
		sz = [];
		J = [];
	end
	nc = J;
end

%% function #3
% % % modified 4/11/2013 both end and 1st
function	[outBinArray, nNeighborsInTimeFreq, nNeighborsTopoSpatial] = ...
	spatialFilterCluster4Coupling(binArray, neighborDefinition, spfThresh)

% % %   17/09/2014: updated by wp fpr better compatibility
	
	% % % prepare para
	sz = size(binArray);	%expect(f, t, node, node)
	nNodes = length(neighborDefinition);
	
	% % % numbers of possible spatial neighbors	
	nNeighborsPossiblePerNode = zeros(nNodes, 1);
	for k = 1 : nNodes
		nNeighborsPossiblePerNode(k) = length(neighborDefinition{k}) - 1;	%spatial neighbor per node
	end
	nNeighborsPossiblePerEdge = bsxfun(@plus, nNeighborsPossiblePerNode, nNeighborsPossiblePerNode');
	nNeighborsPossiblePerEdge = permute(nNeighborsPossiblePerEdge, [3 4 1 2]);
	
	% % % numbers of possible time-frequency neighbors	
	nNeighborsPossiblePerTimeFreq = zeros(sz(1 : 2)) + 4;
% 	nNeighborsPossiblePerTimeFreq([1, end], :) = nNeighborsPossiblePerTimeFreq([1, end], :) - 1;
% 	nNeighborsPossiblePerTimeFreq(:, [1, end]) = nNeighborsPossiblePerTimeFreq(:, [1, end]) - 1;
	nNeighborsPossiblePerTimeFreq(1, :) = nNeighborsPossiblePerTimeFreq(1, :) - 1;
	nNeighborsPossiblePerTimeFreq(:, 1) = nNeighborsPossiblePerTimeFreq(:, 1) - 1;
	nNeighborsPossiblePerTimeFreq(end, :) = nNeighborsPossiblePerTimeFreq(end, :) - 1;
	nNeighborsPossiblePerTimeFreq(:, end) = nNeighborsPossiblePerTimeFreq(:, end) - 1;

	% % % real number of Neighbors	
	nNeighborsInTimeFreq = tfFilter(binArray);
	nNeighborsTopoSpatial = topoFilter(binArray, neighborDefinition);
	
	% % % sum together and apply filter	
	nNeighborsPossibleAll = bsxfun(@plus, nNeighborsPossiblePerEdge, nNeighborsPossiblePerTimeFreq);
	outBinArray = binArray & ((nNeighborsInTimeFreq + nNeighborsTopoSpatial) > spfThresh * nNeighborsPossibleAll);
end

%% function #4
function	nNeighborsInTimeFreq = tfFilter(binArray)

	sz = size(binArray);	%expect ([f, t, nodes, nodes])
	tmpArray = zeros([sz, 4]);
	% % % shift array to 	
	tmpArray(2 : end, :, :, :, 1) = binArray(1 : end - 1, :, : ,:);
	tmpArray(1 : end - 1, :, :, :, 2) = binArray(2 : end, :, : ,:);
	tmpArray(:, 2 : end, :, :, 3) = binArray(:, 1 : end - 1, : ,:);
	tmpArray(:, 1 : end - 1, :, :, 4) = binArray(:, 2 : end, : ,:);
	nNeighborsInTimeFreq = sum(tmpArray, 5);
	% % % set self connection as no neigbors	
	for ig = 1 : sz(4)
		nNeighborsInTimeFreq(:, :, ig, ig) = 0;
	end
end

%% function #5
function	nNeighborsTopoSpatial = topoFilter(binArray, neighborDefinition)

	sz = size(binArray);	%expect(f, t, node, node)
	nNeighborsTopoSpatial = zeros(sz);
	if sz(4) ~= length(neighborDefinition)
		error('inconsistent grid number in connection and neighbor definition!');
	end
	
	for iNode = 1 : sz(4)
		idx = neighborDefinition{iNode};	%all neighbors
		idx(idx == iNode) = [];	%exclude itself
		nNeighborsTopoSpatial(:, :, :, iNode) = nNeighborsTopoSpatial(:, :, :, iNode) + sum(binArray(:, :, :, idx), 4);
		nNeighborsTopoSpatial(:, :, iNode, :) = nNeighborsTopoSpatial(:, :, iNode, :) + sum(binArray(:, :, idx, :), 3);
		nNeighborsTopoSpatial(:, :, iNode, iNode) = 0;	%auto connection is set as 0
	end
% 	nNeighborsTopoSpatial(nNeighborsTopoSpatial < 0) = 0;
end

%% function #6
function [c, nc] = findCluster4CouplingNode(con)

% % %   17/09/2014: updated by wp fpr better compatibility
	

	sz =size(con);
	for ch1 = 1 : sz(4) - 1
		con(:, :, ch1, ch1 : sz(3)) = 0;
	end
	[x1, x2, x3, x4] = ind2sub(sz, find(con));
	nds = length(x1);
% 	fprintf('Total number of rooms :%d of %d.\n', nds, numel(con));	%###
	fg = true(nds, 1);
	ids = [x3, x4, x1, x2];
	nc = 0;
	c = [];
% 	fprintf('3rd floor!\n'); %###
	for ig = 1 : nds
% 		fprintf('room #%d!\n', ig); %###
		if fg(ig)
			cid = [];
			[cid, fg] = searchClusterLoopSparse(ids, ig, cid, fg);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :);
				c{nc}.nd = unique(ids(cid, 1:2));
				c{nc}.nn = length(c{nc}.nd);
				c{nc}.id = sub2ind(sz, x1(cid), x2(cid), x3(cid), x4(cid));
			end
		end
	end
% 	fprintf('4th floor!\n'); %###
end

%% function #7
function [y, flag] = searchClusterLoopSparse(c4, x, y, flag)

	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			s1 = find(sum(bsxfun(@eq, c4, c4(x(k), :)), 2) >= 3);
			s2 = find(sum(bsxfun(@eq, c4, c4(x(k), [2 1 3 4])), 2) >= 3);
			s = unique([s1; s2]);
			s(s==x(k)) = [];
			d = sum(abs(bsxfun(@minus, c4(s, 3:4), c4(x(k), 3:4))), 2) > 1;
			s(d) = [];
			x2 = [x2; s];
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparse(c4, x2, y, flag);
	end
end

%% function #8
function [c, nc] = findCluster4CouplingEdge(con, nb)


% % %   17/09/2014: updated by wp fpr better compatibility
	
	sz =size(con);
	for ch1 = 1 : sz(4) - 1
		con(:, :, ch1, ch1 : sz(3)) = 0;
	end
	[x1, x2, x3, x4] = ind2sub(sz, find(con));
	nds = length(x1);
% 	fprintf('Total number of rooms :%d of %d.\n', nds, numel(con));	%###
	fg = true(nds, 1);
	ids = [x3, x4, x1, x2];
	nc = 0;
	c = [];
% 	fprintf('3rd floor!\n'); %###
	for ig = 1 : nds
% 		fprintf('room #%d!\n', ig); %###
		if fg(ig)
			cid = [];
			[cid, fg] = searchClusterLoopSparseNB(ids, ig, cid, fg, nb);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :);
				c{nc}.nd = unique(ids(cid, 1:2));
				c{nc}.nn = length(c{nc}.nd);
				c{nc}.id = sub2ind(sz, ids(cid, 3), ids(cid, 4), ids(cid, 1), ids(cid, 2));
			end
		end
	end
% 	fprintf('4th floor!\n'); %###
end

%% function #9
function [y, flag] = searchClusterLoopSparseNB(c4, x, y, flag, nb)

	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			s1 = find(sum(bsxfun(@eq, c4(:, 3:4), c4(x(k), 3:4)), 2) == 2);
			if ~isempty(s1)
				for k1 = 1 : 2
					for k2 = 1 : 2
						s2 = find(ismember(c4(s1, k2), c4(x(k), k1)));
						if ~isempty(s2)
							s3 = find(ismember(c4(s1(s2), 3 - k2), nb{c4(x(k), 3 - k1)}));
							if ~isempty(s3)
								x2 = [x2; s1(s2(s3))];
							end
						end
					end
				end
			end
			clear s1 s2 s3 k1 k2;
			for k1 = 1 : 2
				id = [1, 2, k1 + 2];
				s1 = find(sum(bsxfun(@eq, c4(:, id), c4(x(k), id)), 2) == 3);
				if ~isempty(s1)
					s2 = abs(bsxfun(@minus, c4(s1, 5 - k1), c4(x(k), 5 - k1))) == 1;
					if ~isempty(s2)
						x2 = [x2; s1(s2)];
					end
				end
			end
			clear s1 s2 k1;
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparseNB(c4, x2, y, flag, nb);
	end
end

%% function #10
function [th, H] = distr2thresh4H0 (H0, p)
% % % updated 15/11/2017 by wp, nan was treated as -inf
% % % written 03/09/2014 by wp
% % % H0 shall be Nxm matrix, when N is number of cases and m is numbers
% % % within cases, thus N > m. First column contains the largest values
% % % p is the demanded p value (e.g. 0.05)

	%% 1. check input
	% % % data dimension	
	sz = size(H0);
	if sz(1) < sz(2)
		H0 = H0';
		sz = size(H0);
		fprintf('The input seemed to be in the wrong dimension, changed!\n');
	end
	% % % deal with nan
	H0(isnan(H0)) = -Inf;
	% % % default value
	if nargin < 2 || isempty(p)
		p = 0.05;
	end
	
	%% 2. do it
	% % % split data
	H1 = H0(:, 1);
	H2 = H0(:, 2:sz(2));
	% % % sort	
	[H1, I] = sort(H1, 'descend');
	H2 = H2(I, :);
	% % % first round	
	tmp = round(p * sz(1));
	if tmp && (~isempty(H1))
		th1 = H1(tmp);
		H3 = H2(H2 > th1);
		H = sort([H1; H3], 'descend');
		th = H(round(p * length(H)));
	else
		th = [];
		H = [];
	end
	
end

%% function #11
function tfd = cluster_TF2(ed, nFs, nTs)
	tfd = zeros(nFs, nTs);
	for it = 1 : nTs
		for iq = 1 : nFs
			tfd(iq, it) = length(find(ed(:, 3) == iq & ed(:, 4) == it)); 
		end
	end
end %end of function

%% function #12
function [spd1, spd2] = cluster_SP2(ed, nvxs)
	tmp = ed(:, 1:2);
	spd2 = zeros(nvxs);
	for k =  1 : size(tmp, 1)
		spd2(tmp(k, 1), tmp(k, 2)) = spd2(tmp(k, 1), tmp(k, 2)) + 1;
		spd2(tmp(k, 2), tmp(k, 1)) = spd2(tmp(k, 2), tmp(k, 1)) + 1;
	end
% 	spd2 = spd2 / 2; %- modified 04/01/2016
	spd1 = sum(spd2, 2);
	clear tmp;
end %end of function

%% function #13
function p = cluster_pValue(H, sz)
	[y, I] = min(abs(H - sz));
	p = I ./ length(H);
	clear y I H sz;
end %end of function

%% function #14
function [t2, t1, ttf, te] = ed2tClusterCoupling(ed, t, N, nf, nt)
% % % 07/12/2014	written by wp
% % % 	ed([v1, v2, f, t])
% % %		t(f, t, v, v)

	%% prepare
	if nargin < 5
		nt = max(ed(:, 4));
	end
	if nargin < 4
		nf = max(ed(:, 3));
	end
	if nargin < 3
		N = max(max(ed(:, 1:2)));
	end
	t1 = zeros(N, 1);
	t2 = zeros(N, N);
	n = size(ed, 1);
	te = zeros(n, 1);
	ttf = zeros(nf, nt);
	
	%% work
	for k = 1 : n
		te(k) = t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		t1(ed(k, 1)) = t1(ed(k, 1)) + t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		t1(ed(k, 2)) = t1(ed(k, 2)) + t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		t2(ed(k, 1), ed(k, 2)) = t2(ed(k, 1), ed(k, 2)) + t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		t2(ed(k, 2), ed(k, 1)) = t2(ed(k, 2), ed(k, 1)) + t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		ttf(ed(k, 3), ed(k, 4)) = ttf(ed(k, 3), ed(k, 4)) + ...
			t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
	end
		
end % end of function
