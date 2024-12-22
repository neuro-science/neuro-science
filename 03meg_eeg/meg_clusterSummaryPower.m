%% 9. summary
function [cc, th, H] = meg_clusterSummaryPower(H0, c1, c2, nTs, nFs, nvxs, pth, t)
% % % updated 29/01/2019 move summary part out as stand-alone function
% % % 15/07/2015	updated by wp, t
% % % 26/02/2015	written by wp

	%% prepare
	% % % paras
	if nargin < 8 || isempty(t)
		tflag = false;
	elseif any([size(t, 1) - nFs, size(t, 2) - nTs, size(t, 3) - nvxs])
		str1 = sprintf('%d ', size(t));
		str2 = sprintf(' -vs.- ');
		str3 = sprintf('%d ', [nFs, nTs, nvxs]);
		str4 = sprintf(': wrong size of t!\n');
		str = [str1, str2, str3, str4];
		error(str);
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
	[th, H] = mat_distr2thresh(H0, pth);
	
	% % % select the clusters
	cc = [];
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
	if ct
		sz1 = zeros(ct, 1);
		for k = 1 : ct
			sz1(k) = cc{k}.sz;
		end
		[tmp, I] = sort(abs(sz1), 'descend'); %modeified 09/05/2017
		cc = cc(I);
		clear tmp I sz1;
	end
	
	fprintf('=========Data sorted on %04d-%02d-%02d %02d:%02d:%02d=========\n', round(clock));
end % end of function

%% 10. tf2d with n
function tfd = cluster_TF2(ed, nFs, nTs)
	tfd = zeros(nFs, nTs);
	for it = 1 : nTs
		for iq = 1 : nFs
			tfd(iq, it) = length(find(ed(:, 2) == iq & ed(:, 3) == it)); 
		end
	end
end %end of function

%% 11. TF2d with t
function tft = cluster_TFT2(ed, nFs, nTs, t)
	tft = zeros(nFs, nTs);
	for k = 1 : size(ed, 1)
		tft(ed(k, 2), ed(k, 3)) = tft(ed(k, 2), ed(k, 3)) + t(ed(k, 2), ed(k, 3), ed(k, 1)); 
	end
end %end of function

%% 12. sp 1d with n
function spd1 = cluster_SP1(ed, nvxs)
	tmp = ed(:, 1);
	spd1 = zeros(nvxs, 1);
	for k =  1 : size(tmp, 1)
		spd1(tmp(k)) = spd1(tmp(k)) + 1;
	end
	clear tmp;
end %end of function

%% 13. sp 1d with t
function spt1 = cluster_SPT1(ed, nvxs, t)
	spt1 = zeros(nvxs, 1);
	for k =  1 : size(ed, 1)
		spt1(ed(k, 1)) = spt1(ed(k, 1)) + t(ed(k, 2), ed(k, 3), ed(k, 1));
	end
	clear tmp;
end %end of function

%% 14. cluster p value
function p = cluster_pValue(H, sz)
	[y, I] = min(abs(H - sz));
	p = I ./ length(H);
	clear y I H sz;
end %end of function
