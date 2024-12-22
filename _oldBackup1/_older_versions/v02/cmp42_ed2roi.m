function r = cmp42_ed2roi(ed, rID, nFs, nTs)

% % % written  11/12/14 by wp: edges to rois

	%% prepare
	if nargin < 2
		error('edges and rois are needed!');
	end
	if nargin < 3 || isempty(nFs)
		nFs = max(ed(:, 3));
	end
	if nargin < 4 || isempty(nTs)
		nTs = max(ed(:, 4));
	end

	r.vxs = unique(ed(:, 1:2));
	nv = numel(r.vxs);
	vxs2iv = zeros(max(r.vxs), 1);
	for iv = 1 : nv
		vxs2iv(r.vxs(iv)) = iv;
	end
	
	v = cell(nv, 1);
	r.pr = zeros(nv);
	r.tf3 = zeros(nFs, nTs, nv);
	
	%% check the ed
	for iv = 1 : nv
		id = ed(:, 1) == r.vxs(iv) | ed(:, 2) == r.vxs(iv);
		d = ed(id, :);
		for k = 1 : size(d, 1)
			if d(k, 1) == r.vxs(iv)
				r.pr(vxs2iv(d(k, 2)), iv) = r.pr(vxs2iv(d(k, 2)), iv) + 1;
			else
				r.pr(vxs2iv(d(k, 1)), iv) = r.pr(vxs2iv(d(k, 1)), iv) + 1;
			end
			
			r.tf3(d(k, 3), d(k, 4), iv) = r.tf3(d(k, 3), d(k, 4), iv) + 1;
		end
	end
	r.nvs = length(r.vxs);
	ids = rID(r.vxs);
	rids = unique(ids);
	r.nrs = length(rids);
	r.rvx = cell(r.nrs, 1);
% 	r.rlb = cell(r.nrs, 1);
	r.nps = zeros(r.nrs, 1);
	r.nes = zeros(r.nrs, 1);
	for ir = 1 : r.nrs
		r.rvx{ir} = r.vxs(ids == rids(ir));
% 		r.rlb{ir} = rLb{r.rvx{ir}(1)};
		r.nps(ir) = length(r.rvx{ir});
		r.nes(ir) = sum(sum(r.pr(:, ids == rids(ir))));
	end
	
	[r.nes, I] = sort(r.nes, 'descend');
	r.rvx = r.rvx(I);
% 	r.rlb = r.rlb(I);
	r.nps = r.nps(I);
	
	r.rc = zeros(r.nrs, r.nrs);
	for ir1 = 1 : r.nrs - 1
		[tmp1, s1] = ismember(r.rvx{ir1}, r.vxs);
		r.rc(ir1, ir1) = sum(sum((r.pr(s1, s1))));
		for ir2 = ir1 + 1 : r.nrs
			[tmp2, s2] = ismember(r.rvx{ir2}, r.vxs);
			r.rc(ir1, ir2) = sum(sum(r.pr(s1, s2)));
			r.rc(ir2, ir1) = sum(sum(r.pr(s2, s1)));
		end
	end	
	[tmp1, s1] = ismember(r.rvx{r.nrs}, r.vxs);
	r.rc(r.nrs, r.nrs) = sum(sum((r.pr(s1, s1))));
	
end

