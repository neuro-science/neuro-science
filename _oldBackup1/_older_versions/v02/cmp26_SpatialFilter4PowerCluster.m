function	[outBinArray, nNeighborsInTimeFreq, nNeighborsTopoSpatial] = ...
	cmp26_SpatialFilter4PowerCluster(binArray, neighborDefinition, spfThresh)
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
