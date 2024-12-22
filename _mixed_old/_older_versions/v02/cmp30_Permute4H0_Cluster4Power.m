function [H0, myClusters1, sz1, nc1, myClusters2, sz2, nc2] = ...
	cmp30_Permute4H0_Cluster4Power(data, nb, th1, th2, nClustersPerTest, nTestPerCluster, nblimit, fid)
% % % updated 07/02/2017 to allow bigger subject number (up to 90)
	

	%% data preparison
	% % % para defaults
	if nargin < 2
		error('At least four inputs needed: connection array, neigbor definition.');
	end
	if nargin < 3 || isempty(th1)
		th1 = 2.86;
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
	if nargin < 8 || isempty(fid)
		fid = 1;
	end
	
	
	% % % data [f, t, ch, sb]
	[nFs, nTs, nchs, nSubs] = size(data);
	
	%% initialize
	H0 = zeros(nClustersPerTest, nTestPerCluster);
	if nTestPerCluster > 2^nSubs
		error('The number of subjects is too small!');
	elseif nSubs <= 30
		s0 = uint32(randperm(2.^nSubs));
		s1 = ismember(dec2bin(s0(1 : nTestPerCluster) - 1, nSubs), '1');
		clear s0;
	else
		s1 = rand(nTestPerCluster, nSubs) > 0.5;
	end
	
	%% do computation
	if nargout > 1 %computation only
		con = cmp04_Tvalue(data, 4, 0);
		[myClusters1, sz1, nc1] = cmp29_Cluster4Power (con, nb, th1, th2, nblimit);
		[myClusters2, sz2, nc2] = cmp29_Cluster4Power (-con, nb, th1, th2, nblimit);

	else %statistics
		parfor ts = 1 : nTestPerCluster
% 		for ts = 1 : nTestPerCluster
			tic;

			% % % prepare data		
			cz = bsxfun(@times, data, 2 * (permute(s1(ts, :), [1 3 4 2]) - 0.5));
			con = cmp04_Tvalue(cz, 4, 0);

			% % % do cluster
			[tmp1, sz0, tmp3] = cmp29_Cluster4Power (con, nb, th1, th2);
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
			fprintf(fid, 'Test #%03d done after %7.2f seconds.\n', ts, toc);
		end
	end
	
end %end of function

