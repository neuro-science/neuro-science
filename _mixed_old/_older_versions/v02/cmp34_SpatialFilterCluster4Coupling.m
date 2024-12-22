% % % modified 4/11/2013 both end and 1st
function	[outBinArray, nNeighborsInTimeFreq, nNeighborsTopoSpatial] = ...
	cmp34_SpatialFilterCluster4Coupling(binArray, neighborDefinition, spfThresh)

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
