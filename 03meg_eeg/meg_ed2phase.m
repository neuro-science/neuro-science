function [a, pwr] = meg_ed2phase(ed, tf, trls, delay)
% % % 29/04/2019	updated by wp, add cluster4power, no ABS, version 5
% % % 29/11/2017	updated by wp, time delay added, also trls directed to phase base version 4
% % % 05/04/2017	updated by wp, time delay added, also trls directed to phase base version 4
% % % 06/04/2016	written by wp, version 3
% % % 	ed([v1, v2, f, t])
% % %		tf({t, f}(ch, tr, tp))

	%% prepare
	% % % delays? 
	if nargin < 4 || isempty(delay)
		fprintf('no delays in edge!\n');
	elseif numel(delay) == 1
		fprintf('There is %d time delay in edge!\n', delay);
		ed = bsxfun(@plus, ed, [0 0 0 delay]);
	elseif size(delay) == [1 4]
		fprintf('There is %d delay in edge!\n', delay);
		ed = bsxfun(@plus, ed, delay);
	else
		error('delay was ambiguous!');
	end
	% % % trial selection?
	ntr = size(tf{1, 1}, 2);
	n = size(ed, 1);
	if nargin < 3 || isempty(trls)
		trls = 1 : ntr;
	end
	
	
	%% work
	if size(ed, 2) == 4	%coupling
		% % % initialize
		for k = n : -1 : 1
			d1 = tf{ed(k, 4), ed(k, 3)}(ed(k, 1), :, :);
			d2 = tf{ed(k, 4), ed(k, 3)}(ed(k, 2), :, :);
			tmp1 = abs(d1) .* abs(d2);
			d = (d1 .* conj(d2)) ./ (tmp1);
			d0 = mean(d(:, trls, :), 2);
			a(k, :) = angle(mean(bsxfun(@times, d, conj(d0)), 3));
			tmp2 = mean(tmp1, 3);
			pwr(k, :) = bsxfun(@minus, tmp2, mean(tmp2(:, trls), 2));
			clear tmp1 tmp2;
		end

		if nargout < 2
			pwr = [];
		end
	elseif size(ed, 2) == 3	%power
		% % % initialize
		for k = n : -1 : 1
			d1 = tf{ed(k, 3), ed(k, 2)}(ed(k, 1), :, :);
			d = d1 ./ abs(d1);
			d0 = mean(d(:, trls, :), 2);
			a(k, :) = angle(mean(bsxfun(@times, d, conj(d0)), 3));
		end
		pwr = [];
	else
		error('size of edges are not recognized!');
	end
end % end of function
