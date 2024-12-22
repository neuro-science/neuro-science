function [c, sz] = cmp22_cluster2d (t, th)
% % % written 02/09/2014

	%% 1. get the clusters
	% % % input size
	[nfs, nts] = size(t);
	% % %	thresholding
	tmp = t > th;
	% % % compute
	y = bwconncomp(tmp, 4);
	clear tmp;
	
	% % % 	initialize
	c = cell(y.NumObjects, 1);
	sz = zeros(y.NumObjects, 1);
	
	% % % loop for contents	
	for k = 1 : y.NumObjects
		[c{k}.fid, c{k}.tid] = ind2sub([nfs, nts], y.PixelIdxList{k});
		c{k}.sz = sum(t(y.PixelIdxList{k}));
		sz(k) = c{k}.sz;
	end
	[sz, I] = sort(sz, 'descend');
	c = c(I);
	clear y I k;

end
