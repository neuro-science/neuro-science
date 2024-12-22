function out = cmp50_coupling_TF_Shiftpredictor (data, trls, methods, flag1d, nRep)
% % % rewrote 13/05/14, flag1d added to allow 2d output

	%% 1. check data and flags
	
	% % % check input
	if nargin < 5 || isempty(nRep)
		out.nRep = 1000;	%default 1000 repetation
	else
		out.nRep = nRep;
	end
	
	if nargin < 4
		flag1d = 0;	%default 2D
	end
	
	if nargin < 3
		methods = [];	%default no coherence, only csd
	end
	
	if nargin < 2 || isempty(trls)
		trls = 1 : size(data{1, 1}, 2);	%default all trials
	end
	
	% % % size etc.
	[out.nTs, out.nFs] = size(data);
	out.nChs = size(data{1, 1}, 1);
	out.nTrls = length(trls);
	out.nPairs = out.nChs * (out.nChs - 1) / 2;
	
	% % % sub2ind
	out.id2 = zeros(out.nPairs, 2);
	ct = 0;	% % % counter
	for ch1 = 1 : out.nChs - 1
		for ch2 = ch1 + 1 : out.nChs
			ct = ct + 1;
			out.id2(ct, 1) = ch1;
			out.id2(ct, 2) = ch2;
		end
	end
	out.id1 = sub2ind([out.nChs, out.nChs], out.id2(:, 1), out.id2(:, 2));

	if flag1d
		out.size = [out.nRep, out.nPairs, out.nTs, out.nFs];
	else
		out.size = [out.nRep, out.nChs, out.nChs, out.nTs, out.nFs];
	end
	???Start here!???
	% % % check methods
	mFlag = zeros(3, 1);	%[coh, plv, imc]
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
	for iq = 1 : out.nFs
		for it = 1 : out.nTs
			nTapers = size(data{it, iq}, 3);
			N = out.nTrls * nTapers;
			for ch1 = 1 : out.nChs - 1
				d1 = reshape(data{it, iq}(ch1, trls, :), [out.nChs, N]);
				for ch2 = ch1 + 1 : out.nChs
					for k = 1 : n
					

			% % % pre-calculate
			ps = d .* conj(d);
			p = sum(ps, 2);
			c = d * d';
			if flag1d
				out.csd(:, it, iq) = c(out.id1) / N;
			else
				out.csd(:, :, it, iq) = c / N;
			end
			
			% % % coherence
			if mFlag(1)
				coh = abs(c) ./ sqrt(p * p');
				if flag1d
					out.coh(:, it, iq) = coh(out.id1);
				else
					out.coh(:, :, it, iq) = coh;
				end
			end

			% % % plv		
			if mFlag(2)
				d0 = d ./ sqrt(ps);
				plv = abs(d0 * d0') / N;
				if flag1d
					out.plv(:, it, iq) = plv(out.id1);
				else
					out.plv(:, :, it, iq) = plv;
				end
			end

			% % % imc		
			if mFlag(3)
				imc = imag(c) ./ sqrt(p * p');
				if flag1d
					out.imc(:, it, iq) = imc(out.id1);
				else
					out.imc(:, :, it, iq) = imc;
				end
			end

		end %end of t loop
	end %end of f loop
end %end of function