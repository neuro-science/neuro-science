function r = cmp70_Cluster4PowerInANOVA(data, nb, th1, th2, nClustersPerTest, nTestPerCluster, nblimit, fid, nCores)
% % % written on 08/05/2017, based on the cmp30_Permute4H0_Cluster4Power
	

	%% 1. data preparison
	% % % para defaults
	if nargin < 2
		error('At least four inputs needed: connection array, neigbor definition.');
	end
	if nargin < 3 || isempty(th1)
		th1 = 2.86;	%the threshold F-value for single connection
	end
	if nargin < 4 || isempty(th2)
		th2 = 0.5; %the threshold of spatial filter
	end
	if nargin < 5 || isempty(nClustersPerTest)
		nClustersPerTest = 3; %number clusters within each permutation
	end
	if nargin < 6 || isempty(nTestPerCluster)
		nTestPerCluster = 1000; % number of permutations
	end
	if nargin < 7 || isempty(nblimit)
		nblimit = 1;	%whether the cluster should be contrained by neighbors
	end
	if nargin < 8 || isempty(fid)
		fid = 1;	%file handle for output message (1-screen)
	end
	if nargin < 9 || isempty(fid)
		nCores = 1;	%file handle for output message (1-screen)
	end
	if nargin < 10 || isempty(fname)
		tmp = clock;
		fname = ['~/tmp_cluster_on', num2str(tmp(1), '%04d'), num2str(tmp(2), '%02d'), num2str(tmp(3), '%02d'), num2str(tmp(4), '%02d'), num2str(tmp(5), '%02d'), '_', num2str(round(tmp(6)*1000), '%05d'), '.mat'];	%file handle for output message (1-screen)
	end
	
	% % % data [f, t, ch, sb, lf1, lf2], here lf1 and lf2 are the number of levels for factor 1 and 2
	[nFs, nTs, nvxs, nSubs, lf1, lf2] = size(data);	
	
	%% 2. do computation
	% % % the real data
	[r.F, r.p] = rmANOVA4ThisData(data, 0, nCores);	% % % size: (nFs, nTs, nvxs, 3);
	r.data = cell(3, 1);
	for ic = 1 : 3
		[r.cl{ic}, r.sz{ic}, r.nc{ic}] = cmp29_Cluster4Power (r.F(:, :, :, ic), nb, th1, th2, nblimit);
	end
	save(fname, 'r');
	% % % statistics of null distribution
	r.H0 = zeros(nClustersPerTest, nTestPerCluster, 3);
	for ts = 1 : nTestPerCluster
%	for ts = 1 : nTestPerCluster
		tic;
		% % % prepare data		
		[F, p] = rmANOVA4ThisData(data, 1, nCores);	% % % size: (nFs, nTs, nvxs, 3);

		for ic = 1 : 3
			% % % do cluster
			[tmp1, sz0, tmp3] = cmp29_Cluster4Power (F(:, :, :, ic), nb, th1, th2, nblimit);
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
			r.H0(:, ts, ic) = szn;

			% % % echo message		
			fprintf(fid, 'Test #%03d done after %7.2f seconds.\n', ts, toc);
		end
		save(fname, 'r', '-append');
	end
	
	%% 3. summarize the data
	% % % define some parameters
	PTH = 0.05;
	cstr = {'F1main', 'F2main', 'F1*F2'};
	% % % initialize
	r.th = zeros(3, 1);
	cc = cell(1, 1);
	ct = 0;
	% % % gether the data in a loop
	for ic = 1 : 3
		% % % get the thresh
		[r.th(ic), H] = cmp23_distr2thresh(r.H0(:, :, ic), PTH);
		% % % check the clusters one by one
		for ik = 1 : numel(r.cl{ic})
			if r.cl{ic}{ik}.sz > r.th(ic)
				ct = ct + 1;
				cc{ct} = r.cl{ic}{ik};
				cc{ct}.prefix = cstr{ic};
				cc{ct}.p0 = cluster_pValue(H, r.cl{ic}{ik}.sz);
				if nTs * nFs > 1
					cc{ct}.tfd = cluster_TF2(cc{ct}.ed, nFs, nTs);
					cc{ct}.tff = cluster_TFT2(cc{ct}.ed, nFs, nTs, r.F(:, :, :, ic));
				end
				cc{ct}.spd1 = cluster_SP1(cc{ct}.ed, nvxs);
				cc{ct}.spf1 = cluster_SPT1(cc{ct}.ed, nvxs, r.F(:, :, :, ic));
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
		[tmp, I] = sort(sz1, 'descend');
		r.c = cc(I);
		clear tmp I sz1 cc;
	end
	
	% % % save and notify
	save(fname, 'r', '-append');
	fprintf('=========Data sorted on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));	
end %end of function

function [F, p] = rmANOVA4ThisData(data, rFlag, nCores)
	
	%% 1. check data and reformat
	% % % data [f, t, ch, sb, lf1, lf2], here lf1 and lf2 are the number of levels for factor 1 and 2
	[nFs, nTs, nvxs, nSubs, lf1, lf2] = size(data);	
	% % % reshape the data to meet the requirement
	data = reshape(data, [nFs * nTs * nvxs, nSubs* lf1 * lf2])'; %[num of voxels, num of tests]
	% % % parameters for the case-definition
	sb = repmat(1 : nSubs, [1, lf1 * lf2])'; %subjects
	if rFlag % permute the data for control distribution
		seq = zeros(nSubs, lf1*lf2);
		for k = 1 : nSubs
			seq(k, :) = randperm(lf1*lf2);
		end
	else % the real data
		seq = repmat(1 : lf1*lf2, [nSubs, 1]);
	end
	f1 = mod(seq - 1, lf1) + 1; f1 = f1(:); %factor 1
	f2 = ceil(seq / lf1); f2 = f2(:);	%factor 2
	
	%% 2. compute data
	% % % initialize
	F = zeros(nFs * nTs * nvxs, 3);
	p = zeros(nFs * nTs * nvxs, 3);
	% % % compute: a lot of loops, better use parallel computing
	if nCores > 1
		parfor iv = 1 : nFs * nTs * nvxs
			stats = rm_anova2_modified(data(:, iv), sb, f1, f2, {'F1', 'F2'});
			F(iv, :) = [stats{2, 5}, stats{3, 5}, stats{4, 5}];
			p(iv, :) = [stats{2, 6}, stats{3, 6}, stats{4, 6}];
		end
	else
		for iv = 1 : nFs * nTs * nvxs
			stats = rm_anova2_modified(data(:, iv), sb, f1, f2, {'F1', 'F2'});
			F(iv, :) = [stats{2, 5}, stats{3, 5}, stats{4, 5}];
			p(iv, :) = [stats{2, 6}, stats{3, 6}, stats{4, 6}];
		end
	end
	
	%% 3. output: reshape the results
	F = reshape(F, nFs, nTs, nvxs, 3);
	p = reshape(p, nFs, nTs, nvxs, 3);
end


function tfd = cluster_TF2(ed, nFs, nTs)
	tfd = zeros(nFs, nTs);
	for it = 1 : nTs
		for iq = 1 : nFs
			tfd(iq, it) = length(find(ed(:, 2) == iq & ed(:, 3) == it)); 
		end
	end
end %end of function

function tff = cluster_TFT2(ed, nFs, nTs, F)
	tff = zeros(nFs, nTs);
	for k = 1 : size(ed, 1)
		tff(ed(k, 2), ed(k, 3)) = tff(ed(k, 2), ed(k, 3)) + F(ed(k, 2), ed(k, 3), ed(k, 1)); 
	end
end %end of function

function spd1 = cluster_SP1(ed, nvxs)
	tmp = ed(:, 1);
	spd1 = zeros(nvxs, 1);
	for k =  1 : size(tmp, 1)
		spd1(tmp(k)) = spd1(tmp(k)) + 1;
	end
	clear tmp;
end %end of function

function spf1 = cluster_SPT1(ed, nvxs, F)
	spf1 = zeros(nvxs, 1);
	for k =  1 : size(ed, 1)
		spf1(ed(k, 1)) = spf1(ed(k, 1)) + F(ed(k, 2), ed(k, 3), ed(k, 1));
	end
	clear tmp;
end %end of function

function p = cluster_pValue(H, sz)
	[y, I] = min(abs(H - sz));
	p = I ./ length(H);
	clear y I H sz;
end %end of function