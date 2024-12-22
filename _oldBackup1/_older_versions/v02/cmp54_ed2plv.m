function [plva, plza] = cmp54_ed2plv(ed, tf, trls)
% % % 06/04/2016	written by wp, version 3
% % % 	ed([v1, v2, f, t])
% % %		tf({t, f}(ch, tr, tp))

	%% prepare
	if nargin < 3
		trls = size(tf{1, 1}, 2);
	end
	nTrials = length(trls);
	n = size(ed, 1);

	plva = zeros(n, 1);
	plza = zeros(n, 1);

	%% work
	for k = 1 : n
		d1 = tf{ed(k, 4), ed(k, 3)}(ed(k, 1), trls, :);
		d2 = tf{ed(k, 4), ed(k, 3)}(ed(k, 2), trls, :);
		csd = mean(d1 .* conj(d2), 3); %[1grd, trl]
		pw1 = mean(d1 .* conj(d1), 3);	%power [1grd, tr]
		pw2 = mean(d2 .* conj(d2), 3);	%power [1grd, tr]
		df = n * size(d1, 3) - 1;	%
		
		plv0 = csd ./ sqrt(pw1 .* pw2); %[1grd, trl]
		s0 = sum(plv0, 2); %1
		plva(k) = abs(s0 / nTrials); %1
		plza(k) = 23 / 20 * (sqrt(log(1 - plva(k).^2) * (2 - df)) - 23 / 20); %1
	end

end % end of function
