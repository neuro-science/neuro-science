%% function #1
function cc = meg_cluster_pwr_merge(c1, c2, dbFlag)
	% % % written 28/01/2018 by wp : merge smaller clusters together
	% % % check sign	
	if c1.sgn ~= c2.sgn
		fprintf('Cannot be merged: different signs!\n');
		cc = [];
		return;
	else
		cc.sgn = c1.sgn;
	end
	if nargin < 3 || isempty(dbFlag)
		dbFlag = false;
	end
	% % % data merge-able
	cc.ne = c1.ne + c2.ne;
	cc.nn = c1.nn + c2.nn;
	cc.ed = cat(1, c1.ed, c2.ed);
	cc.nd = cat(1, c1.nd, c2.nd);
	cc.sz = c1.sz + c2.sz;
	cc.tfd = c1.tfd + c2.tfd;
	cc.tft = c1.tft + c2.tft;
	cc.spd1 = c1.spd1 + c2.spd1;
	cc.spt1 = c1.spt1 + c2.spt1;
	% % % data not merge-able
	if ~isfield(c1, 'merge') && ~isfield(c2, 'merge')
		cc.merge.p0 = [c1.p0; c2.p0];
		cc.merge.tag = {c1.tag; c2.tag};
		cc.merge.ct = 2;
	elseif ~isfield(c2, 'merge') && isfield(c1.merge, 'ct') && ~isempty(c1.merge.ct)
		cc.merge = c1.merge;
		cc.merge.ct = c1.merge.ct + 1;
		cc.merge.tag{cc.merge.ct} = c2.tag;
		cc.merge.p0(cc.merge.ct) = c2.p0;
	elseif ~isfield(c1, 'merge') && isfield(c2.merge, 'ct') && ~isempty(c2.merge.ct)
		cc.merge = c2.merge;
		cc.merge.ct = c2.merge.ct + 1;
		cc.merge.p0 = [c1.p0; c2.merge.p0];
		cc.merge.tag = cat(1, c1.tag, c2.merge.tag);
	elseif isfield(c1.merge, 'ct') && ~isempty(c1.merge.ct) && isfield(c2.merge, 'ct') && ~isempty(c2.merge.ct) 
		cc.merge.ct = c1.merge.ct + c2.merge.ct;
		cc.merge.p0 = [c1.merge.p0; c2.merge.p0];
		cc.merge.tag = cat(1, c1.merge.tag, c2.merge.tag);
	else
		fprintf('Merge field odd!\n');
	end
	cc.p0 = min(cc.merge.p0);
	xxx = ~diff(char(cc.merge.tag{:}));
	if size(xxx, 1) > 1
		xxx = all(xxx, 1);
	end
	cc.tag = cc.merge.tag{1}(xxx);
	clear xxx;
	% % % double check, no overlapping nodes
	if dbFlag
		nn = numel(unique(cc.nd));
		ne = size(unique(cc.ed, 'rows'), 1);
		if (cc.nn - nn) || (cc.ne - ne)
			fprintf('\n\nError in merged size!!!\n\n\n');
			cc = [];
			cc{2} = c2;
			cc{1} = c1;
		end
	end
end % end of function
