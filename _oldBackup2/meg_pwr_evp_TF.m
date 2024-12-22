function out = meg_pwr_evp_TF (data, trls, chans)
% % % written 05/10/17, power evenlop correlation

	%% 1. check data and flags
	
	% % % check input
	if nargin < 3 || isempty(chans)
		chans = 1 : size(data{1, 1}, 1);	%default all channels
	elseif islogical(chans)
		chans = find(chans);
	end
	
	if nargin < 2 || isempty(trls)
		trls = 1 : size(data{1, 1}, 2);	%default all trials
	elseif islogical(trls)
		trls = find(trls);
	end
	
	% % % size etc.
	[out.nTs, out.nFs] = size(data);
	out.nChs = length(chans);
	out.nTrls = length(trls);
	out.size = [out.nChs, out.nChs, out.nTs, out.nFs];
	out.data = zeros(out.size) + nan;

	%% 2. compute triple methods at once
	for iq = 1 : out.nFs
		for it = 1 : out.nTs
			% % % prepare data
			d1 = permute(data{it, iq}(chans, trls, :), [1 4 2 3]);
			d2 = permute(data{it, iq}(chans, trls, :), [4 1 2 3]);
			d0 = mean(d2 .* conj(d2), 4);
			d = log(bsxfun(@rdivide,	mean(imag(bsxfun(@times, d1, conj(d2))), 4).^2, d0));
			d0 = log(repmat(d0, [out.nChs 1 1]));
			clear d1 d2;

			d1 = bsxfun(@minus, d, mean(d, 3));
			d2 = bsxfun(@minus, d0, mean(d0, 3));
			out.data(:, :, it, iq) = sum(d1 .* d2, 3) ./ sqrt(sum(d1.^2, 3) .* sum(d2.^2, 3));
		end %end of t loop
	end %end of f loop
end %end of function