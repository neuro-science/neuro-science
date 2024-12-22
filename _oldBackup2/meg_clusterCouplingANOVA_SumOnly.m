%% function #1
function c = meg_clusterCouplingANOVA_SumOnly(cl, H0, clp, F, fname)
	% % % updated 19/02/2019 by wp : adopted from meg_clusterCouplingANOVA
	% % % updated 08/02/2019 by wp : re-write to make it working
	% % % updated 02/07/2018 by wp : small message bug corrected
	% % % updated 14/11/2017 by wp : merge several sub functions
	% % % updated 09/03/2017 by wp : add parallel computation
	% % % updated 26/11/2015 by wp : add output 

	%% 01. check input and set defaults
	% % % check input consistency
	if nargin < 2 || isempty(cl) || isempty(H0)
		error('Inputs needed: clusters(Filename), H0(Filename), parameters[optional]');
		c = [];
		return;
	elseif nargin < 4 && ~ischar(cl)
		error('Inputs needed: clusters(Filename), H0(Filename), parameters[optional]');
		c = [];
		return;
	else
		fprintf('Checking inputs ...\n');
	end
	
	% % % check if inputs 1 is data filename
	if ischar(cl)
		try
			v = load(cl);
			fname = cl;
			cl = v.cl;
			clp = v.clp;
			F = v.F;
		catch
			error('Error in load cluster data!');
		end
	end
		
	% % % check if inputs 2 is H0 filename
	if ischar(H0)
		try
			H0 = load(H0);
			if isvector(H0)
				H0 = permute(reshape(H0, [3 numel(H0)/3]), [2 3 1]);
			end
		catch
			error('Error in load H0 data!');
		end
	end
	
	% % % file name
	if ~isvarname('fname')
		if isfield(clp, 'fname')
			fname = clp.fname;
		elseif isfield(clp, 'fileName')
			fname = clp.fileName;
		elseif isfield(clp, 'filename')
			fname = clp.filename;
		else
			error('Error in checking file name!');
		end
	end

	
	%% 05. summarize the data
	% % % prepare
	cstr = {'F1M', 'F2M', 'F1*F2'};
	c = cell(2, 1); %assume not empty
	ct = 0;
	for ic = 3 : -1 : 1
		[th(ic), H] = distr2thresh4H0(H0(:, :, ic), clp.th3);
		for ik = 1 : numel(cl{ic})
			if cl{ic}{ik}.sz > th(ic)
				ct = ct + 1;
				c{ct} = cl{ic}{ik};
				c{ct}.tag = [cstr{ic}, num2str(ik, '%02d')];
				c{ct}.p0 = cluster_pValue(H, cl{ic}{ik}.sz);
				c{ct}.tfd = cluster_TF2(c{ct}.ed, clp.nFs, clp.nTs);
				[c{ct}.spd1, c{ct}.spd2] = cluster_SP2(c{ct}.ed, clp.nVs);
				[c{ct}.spf2, c{ct}.spf1, c{ct}.tff, c{ct}.fe] = ed2tClusterCoupling(c{ct}.ed, F(:, :, :, :, ic), clp.nVs, clp.nFs, clp.nTs);
			end
		end
	end
	
	% % %  sort the clusters
	if ct
		for k = ct : -1 : 1
			sz1(k) = c{k}.sz;
		end
		[tmp, I] = sort(sz1, 'descend');
		c = c(I);
		clear tmp I sz1 k;
	else
		c = [];
	end
	fprintf('=========Data sorted on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));
	
	% % % save the data
	save(fname, 'c', 'th', 'H0', '-append');
	
end %end of function

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
