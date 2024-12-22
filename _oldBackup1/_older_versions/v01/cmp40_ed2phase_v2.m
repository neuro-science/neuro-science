function a = cmp40_ed2phase2(ed, tf, trls)
% % % 08/12/2014	written by wp, version2
% % % 	ed([v1, v2, f, t])
% % %		tf({t, f}(ch, tr, tp))

	%% prepare
	n = size(ed, 1);
	ntr = size(tf{1, 1}, 2);
	a = zeros(n, ntr);
	
	%% work
	for k = 1 : n
		d1 = tf{ed(k, 4), ed(k, 3)}(ed(k, 1), :, :);
		d2 = tf{ed(k, 4), ed(k, 3)}(ed(k, 2), :, :);
		d = (d1 .* conj(d2)) ./ (abs(d1) .* abs(d2));
		a0 = angle(mean(d(1, trls, :), 2));
		ax = bsxfun(@minus, angle(d), a0);
		ay = mod(ax, 2 * pi) ./ (2 * pi);
		a(k, :) = mean(1 - abs(2 * (ay - 0.5)), 3);
	end
	
	
	
end % end of function
