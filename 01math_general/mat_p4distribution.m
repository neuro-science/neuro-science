%% function #13
function p = mat_p4distribution(d, d0)
	d0 = sort(d0, 'descend');
	[y, I] = min(abs(bsxfun(@minus, d0, d)));
	p = I ./ size(d0, 1);
	clear y I H sz;
end %end of function
