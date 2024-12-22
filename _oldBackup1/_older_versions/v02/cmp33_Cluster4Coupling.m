function [myClusters, sz, nc] = cmp33_Cluster4Coupling (conArray, neighborDefinition, conThresh, spfThresh, method)

% % %   17/09/2014: updated by wp fpr better compatibility

% xx_coherence_clustering will find clusters based on the 
% connection array, which defines cnnection strength between locations
% in various time and frequency. The algirthm came from the 2011 Neuron
% paper by Hipp, J.F. Engel, A.K. and Siegel, M.

% Use as
%   myClusters = xx_coherence_clustering (conArray, neighborDefinition, conThresh, spfThresh)
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
%   08/2013: First Formal Draft
%   10/2013: improved Draft to allow more strategies
%   05/2014: improved Draft for more compact sub function

	
	%% preparision
	if nargin < 2
		error('At least two inputs needed: connection array, neigbor definition and two thresholds.');
	elseif nargin < 3
		spfThresh = 0.5;
	elseif nargin < 4
		conThresh = 2.86;	% 
	elseif nargin < 5
		method = 'node';
	end
	
	%% thresholding to binary
	binArray = conArray > conThresh;
	
	%% spatial filtering
	binArray = 	cmp34_SpatialFilterCluster4Coupling (binArray, neighborDefinition, spfThresh);
	%% find cluster
	switch lower(method)
		case 'node'
			[myClusters, nc] = cmp35_findCluster4CouplingNode(binArray);
		case 'edge'
			[myClusters, nc] = cmp36_findCluster4CouplingEdge(binArray, neighborDefinition);
		otherwise
			error('Current supported methods are <node> and <edge> only!');
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
		J = find(nn > 10);
		clear I tmp;
	else
		sz = [];
		J = [];
	end
	nc = J;
end