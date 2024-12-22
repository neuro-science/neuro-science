function [out0, out1] = cmp51_coupling_TF_Shift1predictor (data, trls, methods, flag1d)
% % % written by wp 09/12/15, add shift predictor as controls

	%% 1. check data and flags
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
	[nTs, nFs] = size(data);
	nChs = size(data{1, 1}, 1);
	nTrls = length(trls);
	data1 = cell(size(data));
	
	%% 2. check whether do-able and transform data if yes
	if nChs > nTrls
		error('Number of channels more than number of trials, not supported yet!');
	else
		for iq = 1 : nFs
			for it = 1 : nTs
				nTapers = size(data{it, iq}, 3);
				data1{it, iq} = zeros(nChs, nTrls, nTapers);
				for ch = 1 : nChs
% 					trls1 = trls([ch : end, 1 : ch - 1]);
					trls1 = trls(randperm(nTrls)); %%% this can be improved
					data1{it, iq}(ch, :, :) = data{it, iq}(ch, trls1, :);
				end
			end
		end
	end
	
	%% 3. do it
	out0 = cmp25_coupling_TF(data, trls, methods, flag1d);
	out1 = cmp25_coupling_TF(data1, [], methods, flag1d);
	out1.rd_data = data1;
	
end %end of function

