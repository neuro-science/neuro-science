function [cc, th, H] = cmp44_ClusterSummary(H0, c1, c2, nTs, nFs, nvxs, pth)
% % % 28/03/2017	updated by wp, sort cluster with bigger in the biginning.
% % % 04/01/2016	updated by wp, spd2 would not be divided by 2.
% % % 26/02/2015	written by wp

	%% prepare
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
		nTs = 17;
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
				[cc{ct}.spd1, cc{ct}.spd2] = cluster_SP2(cc{ct}.ed, nvxs);
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
		cc = cc(I);
		clear tmp I sz1;
	end
	
	fprintf('=========Data sorted on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));
end % end of function

function tfd = cluster_TF2(ed, nFs, nTs)
	tfd = zeros(nFs, nTs);
	for it = 1 : nTs
		for iq = 1 : nFs
			tfd(iq, it) = length(find(ed(:, 3) == iq & ed(:, 4) == it)); 
		end
	end
end %end of function

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

function p = cluster_pValue(H, sz)
	[y, I] = min(abs(H - sz));
	p = I ./ length(H);
	clear y I H sz;
end %end of function
