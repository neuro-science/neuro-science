function out = meg_coupling_TF_KtoN (dataK, dataN, trls, methods)
% % % written 21/02/19, K2N data

	%% 1. check data and flags	
	if nargin < 4
		methods = [];	%default no coherence, only csd
	end
	
	if nargin < 3 || isempty(trls)
		trls = 1 : size(dataN{1, 1}, 2);	%default all trials
	else
		trls = trls(:); %mod 2018-12-12, convert multi dimension to vector
	end
	
	% % % size etc.
	[out.nTs, out.nFs] = size(dataN);
	tmp = size(dataK);
	if any([out.nTs ~= tmp(1), out.nFs ~= tmp(2)])
		fprintf('data size do not match: check number of time and frequencies!\n');
		return;
	end
	out.K = size(dataK{1, 1}, 1);
	out.N = size(dataN{1, 1}, 1);
	out.nTrls = numel(trls);
	out.size = [out.K, out.N, out.nTs, out.nFs];
	
	% % % check methods
	mFlag = false(3, 1);	%[coh, plv, imc]
	if any(strcmp('coh', methods))
		mFlag(1) = 1;
		fprintf('Coherence will be computed!\n');
		out.coh = zeros(out.size);
	end
	if any(strcmp('plv', methods))
		mFlag(2) = 1;
		fprintf('Phase locking value will be computed!\n');
		out.plv = zeros(out.size);
	end
	if any(strcmp('imc', methods))
		mFlag(3) = 1;
		out.imc = zeros(out.size);
		fprintf('Imaginary coherence will be computed!\n');
	end
	out.csd = zeros(out.size);

	%% 2. compute triple methods at once
	fprintf('Computing...\n');
	for iq = 1 : out.nFs
		for it = 1 : out.nTs
			% % % prepare data
			nTapers = size(dataN{it, iq}, 3);
			tmp_num = out.nTrls * nTapers;
			dN = reshape(dataN{it, iq}(:, trls, :), [out.N, tmp_num]);
			dK = reshape(dataK{it, iq}(:, trls, :), [out.K, tmp_num]);
			pK = dK .* conj(dK);
			pN = dN .* conj(dN);
			spKN = sqrt(sum(pK, 2) * sum(pN, 2)');
			c = dK * dN';
			out.csd(:, :, it, iq) = c / tmp_num;

			if mFlag(1)
				out.coh(:, :, it, iq) = abs(c) ./ spKN;
				fprintf('c.');
			end
			if mFlag(2)
				d1 = dK ./ sqrt(pK);
				d2 = dN ./ sqrt(pN);
				out.plv(:, :, it, iq) = abs(d1 * d2') / tmp_num;
				fprintf('p.');
			end
			if mFlag(3)
				out.imc(:, :, it, iq) = imag(c) ./ spKN;
				fprintf('i.');
			end
		end %end of t loop
	end %end of f loop
	fprintf('\nComputed!\n');
	
end %end of function