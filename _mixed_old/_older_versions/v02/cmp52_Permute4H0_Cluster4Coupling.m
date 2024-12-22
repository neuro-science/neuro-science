function [H0, myClusters1, sz1, nc1] = ...
	cmp52_Permute4H0_Cluster4Coupling(data, nb, th1, th2, nClustersPerTest, nTestPerCluster, method, oid)
	

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
	if nargin < 7 || isempty(method)
		method = 'node';
	end
	if nargin < 8 || isempty(oid)
		oid = 1;
	end
	
	
	% % % data [f, t, ch, sb]
	[nFs, nTs, nchs, nchs, nSubs] = size(data);
	
	%% do data
	if nargout > 1
		tic;
		con = cmp04_Tvalue(data, 5, 0);
		
		% % % Bug Test
		% 		save('~/3RL/3grp/tmp_c.mat', 'con');
		% % % Bug Test

		fprintf(oid, '\n%s =========Data computed in on %04d-%02d-%02d %02d:%02d:%02d=========\n', method, round(clock));
		[myClusters1, sz1, nc1] = cmp33_Cluster4Coupling (con, nb, th1, th2, method);
		fprintf(oid, '\n%s =========Data computed out on %04d-%02d-%02d %02d:%02d:%02d=========\n', method, round(clock));
		toc;
		H0 = [];
	else
	
		%% initialize
		H0 = zeros(nClustersPerTest, nTestPerCluster);
		s0 = randperm(2.^nSubs);
		s1 = ismember(dec2bin(s0(1 : nTestPerCluster), nSubs), '1');
		clear s0;

		%% do null in loop
		for ts = 1 : nTestPerCluster
			tic;

			% % % prepare data		
			cz = bsxfun(@times, data, 2 * (permute(s1(ts, :), [1 3 4 5 2]) - 0.5));
			con = cmp04_Tvalue(cz, 5, 0);
			cz = [];

			% % % do cluster
			[tmp1, sz0, tmp3] = cmp33_Cluster4Coupling (con, nb, th1, th2, method);
			tmp1 = [];
			tmp2 = [];

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
			try 
				fprintf(oid, 'Test #%03d done after %7.2f seconds.\n', ts, toc);
			catch ME
				fprintf([ME.message, 'Test #%03d done after %7.2f seconds.\n'], ts, toc);
			end			
		end
	end
	
end %end of function

