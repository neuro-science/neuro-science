function [myClusters, sz, nc] = cmp29_Cluster4Power (conArray, neighborDefinition, conThresh, spfThresh, nblimit)
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
	binArray = cmp26_SpatialFilter4PowerCluster(binArray, neighborDefinition, spfThresh);
	
	%% find cluster
	if nblimit
		[myClusters, nc] = cmp27_findCluster4PowerNeighborRestrained (binArray, neighborDefinition);
	else
		[myClusters, nc] = cmp28_findCluster4PowerNeighborFree (binArray);
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