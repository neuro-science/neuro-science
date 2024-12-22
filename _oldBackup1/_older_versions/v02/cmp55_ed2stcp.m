function [stcp, stcp0] = cmp55_ed2stcp(ed, tf, trls)
% % % 06/04/2016	written by wp, version 3
% % % 	ed([v1, v2, f, t])
% % %		tf({t, f}(ch, tr, tp))

	%% prepare
	if nargin < 3
		trls = 1 : size(tf{1, 1}, 2);
	end
	nTrials = length(trls);
	n = size(ed, 1);
	stcp = zeros(n, nTrials);
	stcp0 = zeros(n, nTrials);
	
	%% work
	for k = 1 : n
		d1 = tf{ed(k, 4), ed(k, 3)}(ed(k, 1), trls, :);
		d2 = tf{ed(k, 4), ed(k, 3)}(ed(k, 2), trls, :);
		csd = mean(d1 .* conj(d2), 3);	%[1grd, trl]
		pw1 = mean(d1 .* conj(d1), 3);	%power [1grd, tr]
		pw2 = mean(d2 .* conj(d2), 3);	%power [1grd, tr]
		df = nTrials * size(d1, 3) - 1;	%
		
		plv0 = csd ./ sqrt(pw1 .* pw2); %[1grd, trl]
		s0 = sum(plv0, 2); %1
		plva = abs(s0 / nTrials); %1
		plza = 23 / 20 * (sqrt(log(1 - plva.^2) * (2 - df)) - 23 / 20); %1

		plvi = abs(bsxfun(@minus, s0, plv0) / (nTrials - 1)); %[1, trl]
		plzi = 23 / 20 * (sqrt(log(1 - plvi.^2) * (3 - df)) - 23 / 20); %[1, trl]

		stcp(k, :) = bsxfun(@minus, nTrials * plza, (nTrials - 1) * plzi);
		stcp0(k, :) = bsxfun(@minus, nTrials * plva, (nTrials - 1) * plvi);
	end

end % end of function
