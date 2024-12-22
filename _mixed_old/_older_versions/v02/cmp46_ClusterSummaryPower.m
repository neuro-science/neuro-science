function [cc, th, H] = cmp46_ClusterSummaryPower(H0, c1, c2, nTs, nFs, nvxs, pth, t)
% % % 15/07/2015	updated by wp, t
% % % 26/02/2015	written by wp

	%% prepare
	% % % paras
	if nargin < 8 || isempty(t)
		tflag = false;
	elseif any([size(t, 1) - nvxs, size(t, 2) - nFs, size(t, 3) - nTs])
		fprintf('%d\n', size(t));
		fprintf('%d\n', [nvxs, nFs, nTs]);
		error('wrong size of t!');
	else
		tflag = true;
	end
	
	% % % paras
	if nargin < 7 || isempty(pth)
		pth = 0.05;
	end
	
	if nargin < 6 || isempty(nvxs)
		nvxs = 324;
	end
	
	if nargin < 5 || isempty(nFs)
		nFs = 21;
	end
	
	if nargin < 4 || isempty(nTs)
		nTs = 21;
	end
	
	if nargin < 3
		c2 = [];
	end
	
	if nargin < 2
		error('At leaset 2 inputs are needed!');
	end
	
	% % % data
	c0{1} = c1;
	c0{2} = c2;
	
	%% post-process
	% % % get the thresh
	[th, H] = cmp23_distr2thresh(H0, pth);
	
	% % % select the clusters
	cc = cell(1, 1);
	ct = 0;
	for k = 1 : 2
		for ik = 1 : numel(c0{k})
			if c0{k}{ik}.sz > th
				ct = ct + 1;
				cc{ct} = c0{k}{ik};
				cc{ct}.p0 = cluster_pValue(H, c0{k}{ik}.sz);
				if k > 1.5
					cc{ct}.sgn = -1;
				else
					cc{ct}.sgn = 1;
				end
				cc{ct}.tfd = cluster_TF2(cc{ct}.ed, nFs, nTs);
				cc{ct}.spd1 = cluster_SP1(cc{ct}.ed, nvxs);
				if tflag
					cc{ct}.spt1 = cluster_SPT1(cc{ct}.ed, nvxs, t);
					cc{ct}.tft = cluster_TFT2(cc{ct}.ed, nFs, nTs, t);
				end
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
		[tmp, I] = sort(sz1, 'descend'); %modeified 09/05/2017
		cc = cc(I);
		clear tmp I sz1;
	end
	
	fprintf('=========Data sorted on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));
end % end of function

function tfd = cluster_TF2(ed, nFs, nTs)
	tfd = zeros(nFs, nTs);
	for it = 1 : nTs
		for iq = 1 : nFs
			tfd(iq, it) = length(find(ed(:, 2) == iq & ed(:, 3) == it)); 
		end
	end
end %end of function

function tft = cluster_TFT2(ed, nFs, nTs, t)
	tft = zeros(nFs, nTs);
	for k = 1 : size(ed, 1)
		tft(ed(k, 2), ed(k, 3)) = tft(ed(k, 2), ed(k, 3)) + t(ed(k, 1), ed(k, 2), ed(k, 3)); 
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

function spt1 = cluster_SPT1(ed, nvxs, t)
	spt1 = zeros(nvxs, 1);
	for k =  1 : size(ed, 1)
		spt1(ed(k, 1)) = spt1(ed(k, 1)) + t(ed(k, 1), ed(k, 2), ed(k, 3));
	end
	clear tmp;
end %end of function

function p = cluster_pValue(H, sz)
	[y, I] = min(abs(H - sz));
	p = I ./ length(H);
	clear y I H sz;
end %end of function
