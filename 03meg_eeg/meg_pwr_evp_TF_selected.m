function [cc, ne] = meg_pwr_evp_TF_selected (data, ed, trls, tMove)
% % % written 20/10/17 by wp, power evenlop correlation
% % % ed is supposed to be (v, v, f, t)

	%% 1. check data and flags
	
	% % % check input
	if nargin < 4
		tMove = 0;
	end
	
	if nargin < 3 || isempty(trls)
		trls = 1 : size(data{1, 1}, 2);	%default all trials
	elseif islogical(trls)
		trls = find(trls);
	end
	
	if nargin < 2
		error('You need to provide the data and edges to do computation!');
	end
	
	% % % size etc.
	ne = size(ed, 1);
	cc = zeros(ne, 2);

	%% 2. compute triple methods at once
	for ie = 1 : ne
		% % % prepare data
		d1 = permute(data{ed(ie, 4) + tMove, ed(ie, 3)}(ed(ie, 1), trls, :), [2 3 1]);
		d2 = permute(data{ed(ie, 4) + tMove, ed(ie, 3)}(ed(ie, 2), trls, :), [2 3 1]);
		% % % compute
		y1 = compute_pwr_evp(d1, d2);
		y2 = compute_pwr_evp(d2, d1);
		cc(ie, 1) = (y1 + y2)/2;
		% % % prepare data
		d1 = permute(data{2, ed(ie, 3)}(ed(ie, 1), trls, :), [2 3 1]);
		d2 = permute(data{2, ed(ie, 3)}(ed(ie, 2), trls, :), [2 3 1]);
		% % % compute
		y1 = compute_pwr_evp(d1, d2);
		y2 = compute_pwr_evp(d2, d1);
		cc(ie, 2) = (y1 + y2)/2;
		% % % clean up
		clear d1 d2 y1 y2;
	end
end %end of function


function out = compute_pwr_evp(d1, d2)
	d0 = mean(d2 .* conj(d2), 2);
	d = log(mean(imag(d1 .* conj(d2)), 2).^2 ./ d0);
	d0 = log(d0);
	clear d1 d2;

	d1 = bsxfun(@minus, d, mean(d, 1));
	d2 = bsxfun(@minus, d0, mean(d0, 1));
	out = sum(d1 .* d2, 1) ./ sqrt(sum(d1.^2, 1) .* sum(d2.^2, 1));
end