function [cc, step] = meg_cluster2cluster(spt2, nb, n0, th)
	% % % written 30/05/2018 by wp : cluster to cluster
	M = 20;
	% % % check inputs
	[tmp1, tmp2] = size(spt2);
	if abs(tmp1 - tmp2) > 0.1
		error('connection strength error!');
	else
		n = tmp1;
		clear tmp*;
	end
	
	if abs(numel(nb) - n) > 0.1
		error('neigborhood size error!');
	end
	
	spt1 = sum(spt2, 2);
	S = sum(spt1);
	flag = false(n, 1);
	flag(n0) = true;
	flag(~(spt1)) = true;
	
	N = numel(n0);
	cc = cell(N, 1);
	for k = 1 : N
		cc{k}.id = n0(k);
		cc{k}.cr = 1;
		cc{k}.ds = 0;
		cc{k}.cv = spt2(:, n0(k));
	end
	
	step = 1;
	while sum(spt1(flag))/S < th
		% % % 1st order
		for k = 1 : N
			nb_now = cat(2, nb{cc{k}.id});
			nb_now(flag(nb_now)) = [];
			for kk = 1 : length(nb_now)
				[tmp1, tmp2] = corrcoef(spt2(:, nb_now(kk)), cc{k}.cv);
				tmp(kk) = tmp1(1, 2);
				if tmp(kk) > (1 - step/M)
					cc{k}.id = [cc{k}.id; nb_now(kk)];
					cc{k}.cr = [cc{k}.cr; tmp(kk)];
					cc{k}.ds = [cc{k}.ds; step];
				end
			end
		end
		clear kk;

		for k1 = 1 : N - 1
			for k2 = k1 + 1 : N
				[is, id1, id2] = intersect(cc{k1}.id, cc{k2}.id);
				f1 = [];
				f2 = [];
				if ~isempty(is)
					for kk = 1 : length(is)
						if cc{k1}.cr(id1(kk)) - cc{k1}.ds(id1(kk))/M >= cc{k2}.cr(id2(kk)) - cc{k2}.ds(id2(kk))/M
							f1 = [f1; id1(kk)];
						else
							f2 = [f2; id2(kk)];
						end						
					end
				end
				cc{k2}.id(f2) = [];
				cc{k2}.ds(f2) = [];
				cc{k2}.cr(f2) = [];
				cc{k1}.id(f1) = [];
				cc{k1}.ds(f1) = [];
				cc{k1}.cr(f1) = [];
			end
		end
		for k = 1 : N
			[cc{k}.id, tmp4] = unique(cc{k}.id);
			flag(cc{k}.id) = true;
			cc{k}.ds = cc{k}.ds(tmp4);
			cc{k}.cr = cc{k}.cr(tmp4);
			cc{k}.cv = mean(spt2(:, cc{k}.id), 2);
			cc{k}.sz = sum(spt1(cc{k}.id));
		end
		step = step + 1;
	end
	
		
	
% 	% % % cons
% 	[c1, c2] = corrcoef(spt2);
% % 	c1 = cov(spt2);
% 	c0 = c1 > 0;
% 	
% 	% % % initialize
% 	flag = logical(sum(spt2));
% 	ct = 0;
% 	c = cell(2, 1);
% 	sz = zeros(2, 1);
% 	nc = sz;
% 	
% 	% % % extended neigborhood
% 	rb = cell(n, 1);
% 	for ic = 1 : n
% 		tmp = unique(cat(2, nb{nb{ic}}));
% 		tmp(tmp == ic) = [];
% 		tmp1 = tmp(c0(ic, tmp));
% % 		tmp5 = c1(ic, tmp1) .* sum(spt2(:, tmp1));
% 		tmp5 = c1(ic, tmp1);
% 		[tmp3, tmp2] = sort(tmp5, 'descend');
% 		rb{ic} = tmp1(tmp2);
% 	end
% 	
% 	% % % do search
% 	for ic = 1 : n
% 		if flag(ic)
% 			flag(ic) = false;
% 			ct = ct + 1;
% 			y = ic;
% 			x = rb{ic};
% 			[y, flag] = theRcp(rb, y, x, flag);
% 			c{ct} = sort(y);
% 			sz(ct) = sum(sum(spt2(y, :)));
% 			nc(ct) = numel(y);
% 		end
% 	end
% 	[x1, x2] = sort(abs(sz), 'descend');
% 	sz = sz(x2);
% 	c = c(x2);
% 	nc = nc(x2);
% 	fprintf('\n========\nData was computed @%04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
% 	fprintf('In total %6.2f%% data was explained.\n========\n', sum(sz)*100 / sum(sum(spt2)));
% 	clear x1 x2 y ct ic r flag nb spt2;
	
end % end of function






% function [c, sz, nc, rb] = meg_cluster2cluster(spt2, nb)
% 	% % % written 18/05/2018 by wp : cluster to cluster
% 	
% 	% % % check inputs
% 	[tmp1, tmp2] = size(spt2);
% 	if abs(tmp1 - tmp2) > 0.1
% 		error('connection strength error!');
% 	else
% 		n = tmp1;
% 		clear tmp*;
% 	end
% 	
% 	if abs(numel(nb) - n) > 0.1
% 		error('neigborhood size error!');
% 	end
% 	
% 	% % % cons
% 	[c1, c2] = corrcoef(spt2);
% % 	c1 = cov(spt2);
% 	c0 = c1 > 0;
% 	
% 	% % % initialize
% 	flag = logical(sum(spt2));
% 	ct = 0;
% 	c = cell(2, 1);
% 	sz = zeros(2, 1);
% 	nc = sz;
% 	
% 	% % % extended neigborhood
% 	rb = cell(n, 1);
% 	for ic = 1 : n
% 		tmp = unique(cat(2, nb{nb{ic}}));
% 		tmp(tmp == ic) = [];
% 		tmp1 = tmp(c0(ic, tmp));
% % 		tmp5 = c1(ic, tmp1) .* sum(spt2(:, tmp1));
% 		tmp5 = c1(ic, tmp1);
% 		[tmp3, tmp2] = sort(tmp5, 'descend');
% 		rb{ic} = tmp1(tmp2);
% 	end
% 	
% 	% % % do search
% 	for ic = 1 : n
% 		if flag(ic)
% 			flag(ic) = false;
% 			ct = ct + 1;
% 			y = ic;
% 			x = rb{ic};
% 			[y, flag] = theRcp(rb, y, x, flag);
% 			c{ct} = sort(y);
% 			sz(ct) = sum(sum(spt2(y, :)));
% 			nc(ct) = numel(y);
% 		end
% 	end
% 	[x1, x2] = sort(abs(sz), 'descend');
% 	sz = sz(x2);
% 	c = c(x2);
% 	nc = nc(x2);
% 	fprintf('\n========\nData was computed @%04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
% 	fprintf('In total %6.2f%% data was explained.\n========\n', sum(sz)*100 / sum(sum(spt2)));
% 	clear x1 x2 y ct ic r flag nb spt2;
% 	
% end % end of function
% 
% 
% %% reciprocal 
% function [y, fg] = theRcp(r, y, x, fg)
% 	nn = numel(x);
% 	x2 = [];
% % 						fprintf(' in: ')
% % 						fprintf('y = ')
% % 						fprintf('%4d', y);
% % 						fprintf(';\tx = ');
% % 						fprintf('%4d', x);
% % 						fprintf(';\tx2 = ');
% % 						fprintf('%4d', x2);
% % 						fprintf('.\n');
% 	for k = 1 : nn
% 		if fg(x(k))
% 			for kk = 1 : length(r{x(k)}) 
% 				if ismember(r{x(k)}(kk), y);
% 					fg(x(k)) = false;
% 					y = [y; x(k)];
% 					x2 = unique([x2, r{x(k)}]);
% 					break;
% 				elseif fg(r{x(k)}(kk))
% 					break;
% 				end
% 			end
% 		end
% 	end
% 	y = unique(y);
% 	if ~isempty(x2)
% 		[y, fg] = theRcp(r, y, x2, fg);
% 	end
% % 						fprintf('out: ')
% % 						fprintf('y = ')
% % 						fprintf('%4d', y);
% % 						fprintf(';\tx = ');
% % 						fprintf('%4d', x);
% % 						fprintf(';\tx2 = ');
% % 						fprintf('%4d', x2);
% % 						fprintf('.\n');
% end



% % % function [c, sz, nc] = meg_cluster2cluster(spt2, nb, th)
% % % 	% % % written 18/05/2018 by wp : cluster to cluster
% % % 	
% % % 	% % % check inputs
% % % 	[tmp1, tmp2] = size(spt2);
% % % 	if abs(tmp1 - tmp2) > 0.1
% % % 		error('connection strength error!');
% % % 	else
% % % 		n = tmp1;
% % % 		clear tmp*;
% % % 	end
% % % 	
% % % 	if abs(numel(nb) - n) > 0.1
% % % 		error('neigborhood size error!');
% % % 	end
% % % 	
% % % 	% % % cons
% % % 	[c1, c2] = corrcoef(spt2);
% % % 	c0 = c1 > th;
% % % 	
% % % 
% % % 	flag = logical(sum(c0));
% % % 	ct = 0;
% % % 	c = cell(2, 1);
% % % 	sz = zeros(2, 1);
% % % 	nc = sz;
% % % 	
% % % 	for ic = 1 : n
% % % 		if flag(ic)
% % % 			ct = ct + 1;
% % % 			y = [];
% % % 			x = ic;
% % % 			[y, flag] = theRcp(nb, c0, y, x, flag);
% % % 			c{ct} = sort(y);
% % % 			sz(ct) = sum(sum(spt2(y, :)));
% % % 			nc(ct) = numel(y);
% % % 		end
% % % 	end
% % % 	[x1, x2] = sort(abs(sz), 'descend');
% % % 	sz = sz(x2);
% % % 	c = c(x2);
% % % 	fprintf('\n========\nData was computed @%04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
% % % 	fprintf('In total %6.2f%% data was explained.\n========\n', sum(sz)*100 / sum(sum(spt2)));
% % % 	clear x1 x2 y ct ic r rb flag nb spt2;
% % % 	
% % % end % end of function
% % % 
% % % 
% % % %% reciprocal 
% % % function [y, fg] = theRcp(r, c0, y, x, fg)
% % % 	nn = numel(x);
% % % 	x2 = [];
% % % 	for k = 1 : nn
% % % 		if fg(x(k))
% % % 			fg(x(k)) = false;
% % % 			y = [y; x(k)];
% % % 			x2 = r{x(k)}(logical(c0(r{x(k)}, x(k))));
% % % 			if ~isempty(x2)
% % % 				[y, fg] = theRcp(r, c0, y, x2, fg);
% % % 			end
% % % 		end
% % % 	end
% % % end
% % % 


% % % function [c, sz, nc, rCon] = meg_cluster2cluster(spt2, nb, th1, th2)
% % % 	% % % written 14/05/2018 by wp : cluster to cluster
% % % 	
% % % 	% % % check inputs
% % % 	[tmp1, tmp2] = size(spt2);
% % % 	if abs(tmp1 - tmp2) > 0.1
% % % 		error('connection strength error!');
% % % 	else
% % % 		n = tmp1;
% % % 		clear tmp*;
% % % 	end
% % % 	
% % % 	if abs(numel(nb) - n) > 0.1
% % % 		error('neigborhood size error!');
% % % 	end
% % % 	
% % % 	% % % cons
% % % 	isCon = cell(n, 1);
% % % 	for ic = 1 : n
% % % 		isCon{ic} = find(spt2(:, ic));
% % % 	end
% % % 	rCon = zeros(n);
% % % 	
% % % 	% % % compute ratio
% % % 	for ic = 1 : n
% % % 		x = spt2(ic, :);
% % % 		z = isCon{ic};
% % % 		sx = sum(x);
% % % 		for k = 1 : length(nb{ic})
% % % 			xk = spt2(nb{ic}(k), :);
% % % 			zk = isCon{nb{ic}(k)};
% % % 			[tmp1, tmp2, tmp3] = intersect(z, zk);
% % % 			if ~isempty(tmp1)
% % % 				tmp4 = sum(x(z(tmp2))) ./ sx;
% % % 				tmp5 = sum(xk(zk(tmp3))) ./ sum(xk);
% % % 				if rCon(ic, nb{ic}(k)) && rCon(nb{ic}(k), ic)
% % % 					fprintf('\n\n\n========\n%03d vs. %03d: #3\n========\n\n\n', ic, nb{ic}(k))
% % % 				elseif rCon(ic, nb{ic}(k))
% % % % 					fprintf('%03d vs. %03d: #2!\n', ic, nb{ic}(k))
% % % 					rCon(nb{ic}(k), ic) = (tmp4 + tmp5) / 2;
% % % 				elseif rCon(nb{ic}(k), ic)
% % % 					rCon(ic, nb{ic}(k)) = (tmp4 + tmp5) / 2;
% % % % 					fprintf('%03d vs. %03d: #2!\n', ic, nb{ic}(k))
% % % 				else
% % % 					rCon(ic, nb{ic}(k)) = (tmp4 + tmp5) / 2;
% % % % 					fprintf('%03d vs. %03d: #1!\n', ic, nb{ic}(k))
% % % 				end
% % % 			end
% % % 		end
% % % 	end
% % % 	
% % % 	% % % get cluster
% % % 	rc1 = bsxfun(@rdivide, rCon, sum(rCon, 1));
% % % 	rc2 = bsxfun(@rdivide, rCon, sum(rCon, 2));
% % % 	rc1(isnan(rc1)) = 0;
% % % 	rc2(isnan(rc2)) = 0;
% % % 	r0 = rCon > th1;
% % % 	r1 = rc1 > th2;
% % % 	r2 = rc2 > th2;
% % % 	rb = cell(n, 1);
% % % 	r = r0 | r1 | r2;
% % % 	clear r1 r2 r0;
% % % 	flag = logical(sum(r));
% % % 	for ic = 1 : n
% % % 		rb{ic} = find(r(ic, :));
% % % 	end
% % % 	ct = 0;
% % % 	c = cell(2, 1);
% % % 	sz = zeros(2, 1);
% % % 	nc = sz;
% % % 	for ic = 1 : n
% % % 		if flag(ic)
% % % 			flag(ic) = false;
% % % 			ct = ct + 1;
% % % 			y = ic;
% % % 			x = rb{ic};
% % % 			[y, flag] = theRcp(rb, y, x, flag);
% % % 			c{ct} = sort(y);
% % % 			sz(ct) = sum(sum(spt2(y, :)));
% % % 			nc(ct) = numel(y);
% % % 		end
% % % 	end
% % % 	[x1, x2] = sort(abs(sz), 'descend');
% % % 	sz = sz(x2);
% % % 	c = c(x2);
% % % 	fprintf('\n========\nData was computed @%04d-%02d-%02d %02d:%02d:%02d.\n', round(clock));
% % % 	fprintf('In total %6.2f%% data was explained.\n========\n', sum(sz)*100 / sum(sum(spt2)));
% % % 	clear x1 x2 y ct ic r rb flag nb spt2;
% % % 	
% % % end % end of function
% % % 
% % % 
% % % %% reciprocal 
% % % function [y, fg] = theRcp(r, y, x, fg)
% % % 	nn = numel(x);
% % % 	x2 = [];
% % % 	for k = 1 : nn
% % % 		if fg(x(k))
% % % 			fg(x(k)) = false;
% % % 			y = [y; x(k)];
% % % 			x2 = r{x(k)};
% % % 			if ~isempty(x2)
% % % 				[y, fg] = theRcp(r, y, x2, fg);
% % % 			end
% % % 		end
% % % 	end
% % % end
% % % 
	

% % % version yield 1 cluster
% % % function [c, sz] = meg_cluster2cluster(spt2, nb)
% % % 	% % % written 14/05/2018 by wp : cluster to cluster
% % % 	
% % % 	% % % check inputs
% % % 	[tmp1, tmp2] = size(spt2);
% % % 	if abs(tmp1 - tmp2) > 0.1
% % % 		error('connection strength error!');
% % % 	else
% % % 		n = tmp1;
% % % 		clear tmp*;
% % % 	end
% % % 	
% % % 	if abs(numel(nb) - n) > 0.1
% % % 		error('neigborhood size error!');
% % % 	end
% % % 	
% % % 	% % % extend nb
% % % % 	rb = cell(n, 1);
% % % % 	for ic = 1 : n
% % % % 		rb{ic} = unique([nb{ic}, cat(2, nb{nb{ic}})]);
% % % % 		rb{ic}(rb{ic} == ic) = [];
% % % % 	end
% % % 	rb = nb;
% % % 	
% % % 	% % % get cluster
% % % 	flag = logical(sum(spt2));
% % % 	ct = 0;
% % % 	c = cell(2, 1);
% % % % 	z = [];
% % % 	for ic = 1 : n
% % % 		if flag(ic)
% % % 			flag(ic) = false;
% % % 			ct = ct + 1;
% % % 			y = ic;
% % % % 			fprintf('\n%d ', y);
% % % 			x = rb{ic};
% % % 			[y, flag] = theRcp(rb, y, x, flag);
% % % 			c{ct} = sort(y);
% % % 			sz(ct) = sum(sum(spt2(y, :)));
% % % % 			z = [z; y(:)];
% % % % 			fprintf('\n');
% % % % 			fprintf('%d ', y);
% % % % 			fprintf('\n');
% % % 		end
% % % 	end
% % % 	[x1, x2] = sort(abs(sz), 'descend');
% % % 	sz = sz(x2);
% % % 	c = c(x2);
% % % end % end of function
% % % 
% % % 
% % % %% reciprocal 
% % % function [y, fg] = theRcp(r, y, x, fg)
% % % 	nn = numel(x);
% % % 	x2 = [];
% % % 	for k = 1 : nn
% % % 		if fg(x(k))
% % % 			fg(x(k)) = false;
% % % 			y = [y; x(k)];
% % % % 			fprintf('%d ', x(k));
% % % 			x2 = r{x(k)};
% % % 			if ~isempty(x2)
% % % % 				fprintf('\nx2 = ');
% % % % 				fprintf('%d ', x2);
% % % % 				fprintf('\n');
% % % 				[y, fg] = theRcp(r, y, x2, fg);
% % % 			end
% % % 		end
% % % 	end
% % % end

% % % version yield pieces
% % % function [c, rCon] = meg_cluster2cluster(spt2, nb, th)
% % % 	% % % written 14/05/2018 by wp : cluster to cluster
% % % 	
% % % 	% % % check inputs
% % % 	[tmp1, tmp2] = size(spt2);
% % % 	if abs(tmp1 - tmp2) > 0.1
% % % 		error('connection strength error!');
% % % 	else
% % % 		n = tmp1;
% % % 		clear tmp*;
% % % 	end
% % % 	
% % % 	if abs(numel(nb) - n) > 0.1
% % % 		error('neigborhood size error!');
% % % 	end
% % % 	
% % % 	% % % cons
% % % 	isCon = cell(n, 1);
% % % 	for ic = 1 : n
% % % 		tmp = find(spt2(:, ic));
% % % 		isCon{ic} = unique(cat(2, nb{tmp(:)}));
% % % 	end
% % % 	rCon = zeros(n);
% % % 	
% % % 	% % % compute ratio
% % % 	for ic = 1 : n
% % % 		x = spt2(ic, :);
% % % 		z = isCon{ic};
% % % 		sx = sum(x);
% % % 		for k = 1 : length(nb{ic})
% % % 			xk = spt2(nb{ic}(k), :);
% % % 			zk = isCon{nb{ic}(k)};
% % % 			[tmp1, tmp2, tmp3] = intersect(z, zk);
% % % 			if ~isempty(tmp1)
% % % 				tmp4 = sum(x(z(tmp2))) ./ sx;
% % % 				tmp5 = sum(xk(zk(tmp3))) ./ sum(xk);
% % % 				if rCon(ic, nb{ic}(k)) && rCon(nb{ic}(k), ic)
% % % 					fprintf('\n\n\n========\n%03d vs. %03d: #3\n========\n\n\n', ic, nb{ic}(k))
% % % 				elseif rCon(ic, nb{ic}(k))
% % % 					fprintf('%03d vs. %03d: #2!\n', ic, nb{ic}(k))
% % % 					rCon(nb{ic}(k), ic) = (tmp4 + tmp5) / 2 + 1e-3;
% % % 				elseif rCon(nb{ic}(k), ic)
% % % 					rCon(ic, nb{ic}(k)) = (tmp4 + tmp5) / 2 + 1e-3;
% % % 					fprintf('%03d vs. %03d: #2!\n', ic, nb{ic}(k))
% % % 				else
% % % 					rCon(ic, nb{ic}(k)) = (tmp4 + tmp5) / 2 + 1e-3;
% % % 					fprintf('%03d vs. %03d: #1!\n', ic, nb{ic}(k))
% % % 				end
% % % 			end
% % % 		end
% % % 	end
% % % 	fprintf('\n========\nratio was computed @%04d-%02d-%02d %02d:%02d:%02d.\n========\n', round(clock));
% % % 	
% % % 	% % % get cluster
% % % 	r = rCon > th;
% % % 	rb = cell(n, 1);
% % % 	flag = logical(sum(r));
% % % 	for ic = 1 : n
% % % 		rb{ic} = find(r(ic, :));
% % % 	end
% % % 	ct = 0;
% % % 	c = cell(2, 1);
% % % 	for ic = 1 : n
% % % 		if flag(ic)
% % % 			flag(ic) = false;
% % % 			ct = ct + 1;
% % % 			y = ic;
% % % 			x = rb{ic};
% % % 			[y, flag] = theRcp(rb, y, x, flag);
% % % 			c{ct} = y;
% % % 		end
% % % 	end
% % % 	ct
% % % 	size(c)
% % % 	c = c(1 : ct);
% % % 	
% % % end % end of function
% % % 
% % % 
% % % %% reciprocal 
% % % function [y, fg] = theRcp(r, y, x, fg)
% % % 	nn = numel(x);
% % % 	x2 = [];
% % % 	for k = 1 : nn
% % % 		if fg(x(k))
% % % 			fg(x(k)) = false;
% % % 			y = [y; x(k)];
% % % 			x2 = r{x(k)};
% % % 			if ~isempty(x2)
% % % 				[y, fg] = theRcp(r, y, x, fg);
% % % 			end
% % % 		end
% % % 	end
% % % end

