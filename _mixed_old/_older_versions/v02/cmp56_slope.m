function k = cmp56_slope(d, s)
% % % 12/04/2016	written by wp

	sz = size(d);
	d = reshape(d, sz(1), prod(sz(2 : end)));
	k = zeros(size(d)); 
	
	for ik = 1 : sz(1)
		if ik <= s
			id1 = 1;
			di1 = ik - 1;
		else
			id1 = ik - s;
			di1 = s;
		end
		if ik >= sz(1) - s
			id2 = sz(1);
			di2 = sz(1) - ik;
		else
			id2 = ik + s;
			di2 = s;
		end
		k(ik, :) = (d(id2, :) - d(id1, :)) / (di1 + di2 + 1);
	end
	k = reshape(k, sz);
end % end of function
