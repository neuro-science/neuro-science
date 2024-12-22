function s = mat_mnn_cov (x)
	% % % x(trial, chan); s (chan, chan)
	
	sz = size(x);
	x = bsxfun(@minus, x, mean(x, 1));
	c = cov(x);
	pr = mean(diag(c)) * eye(sz(2));
	xx = x.^2;
	p = sum(sum(xx'*xx/(sz(1)-1) - c.^2));
	g = norm(c - pr, 'fro')^2;
	sk = p / g / sz(1);
	sk = max(0, min(1, sk));
	s = sk * pr + (1 - sk) * c;
	clear c sk g p pr x sz xx;
end