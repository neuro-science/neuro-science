function c = cmp40_ed2phase(ed, tf, trls)
% % % 07/12/2014	written by wp
% % % 	ed([v1, v2, f, t])
% % %		tf({t, f}(ch, tr, tp))

	%% prepare
	n = size(ed, 1);
	if nargin < 3
		ntr = size(tf{1, 1}, 2);
		trls = 1 : ntr;
	else
		ntr = length(trls);
	end
	c = zeros(n, ntr);
	%% work
	for k = 1 : n
		d1 = tf{ed(k, 4), ed(k, 3)}(ed(k, 1), trls, :);
		d2 = tf{ed(k, 4), ed(k, 3)}(ed(k, 2), trls, :);
		d = mean(d1 .* conj(d2), 3);
		a0 = angle(mean(d));
		c(k, :) = mod(angle(d) - a0, 2*pi) ./ (2*pi);
	end
	
	
	
end % end of function
