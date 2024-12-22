function pw = cmp41_ed2pwr(ed, tf)
% % % 07/12/2014	written by wp
% % % 	ed([v1, v2, f, t])
% % %		tf({t, f}(ch, tr, tp))

	%% prepare
	n = size(ed, 1);
	ntr = size(tf{1, 1}, 2);
	pw = zeros(n, ntr);
	
	%% work
	for k = 1 : n
		d1 = abs(tf{ed(k, 4), ed(k, 3)}(ed(k, 1), :, :));
		d2 = abs(tf{ed(k, 4), ed(k, 3)}(ed(k, 2), :, :));
		pw(k, :) = mean(d1 .* d2, 3);
	end
	
	
	
end % end of function
