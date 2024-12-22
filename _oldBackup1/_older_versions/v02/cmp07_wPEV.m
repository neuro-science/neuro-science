function w = cmp07_wPEV(data, trlIDs, wFlag)
% % % 07/07/14	written by wp: omega percentage explained variance
% % % 	data: nTimePoints x nTrials
% % %		trlIDs:	trial labels
% % %		wFlag: 1 for correction with MSE
	%% 1. check inputs
	if nargin < 3
		wFlag = 0;
	end
	npts = size(data, 1);
	if iscell(trlIDs)
		ncons = numel(trlIDs);
		nTrls = zeros(ncons, 1);
		for k = 1 : ncons
			nTrls(k) = length(trlIDs{k});
		end
		trl_all = cat(1, trlIDs{:});
	elseif isnumeric(trlIDs)
		sz = size(trlIDs);
		trl_all = trlIDs(:);
		if length(sz) == 2
			ncons = sz(2);
			nTrls = repmat(sz(1), [1, ncons]);
			tmp = cell(ncons, 1);
			for k = 1 : ncons
				tmp{k} = trlIDs(:, k);
			end
			trlIDs = tmp;
			clear tmp;
		else
			error('I cannot understand 3d trial IDs!');
		end
	else
		error('only numeric or cell format allowed for trial IDs!');
	end

	%% 2. compute
	% % % SS_total
	[SS_total, m_total] = my_var(data(:, trl_all));
	if wFlag
		% % % SS_btGrp, MSE
		SS_btGrp	= 0;
		MSE = 0;
		m_con = zeros(npts, ncons);
		for ic = 1 : ncons
			[se1, m_con(:, ic)] = my_var(data(:, trlIDs{ic}));
			MSE = MSE + se1;
			SS_btGrp = SS_btGrp + nTrls(ic) .* (m_con(:, ic) - m_total).^2;
		end
		% % % output
		w = (SS_btGrp - (ncons - 1) * MSE) ./ (SS_total + MSE);
% 		w1 = w * x;
	else
		% % % SS_btGrp, MSE
		SS_btGrp	= 0;
		m_con = zeros(npts, ncons);
		for ic = 1 : ncons
			m_con(:, ic) = mean(data(:, trlIDs{ic}), 2);
			SS_btGrp = SS_btGrp + nTrls(ic) .* (m_con(:, ic) - m_total).^2;
		end
		% % % output
		w = SS_btGrp ./ SS_total;
	end		
end %end of function

function [v, m] = my_var(data)	%var for 2d data
	m = mean(data, 2);
	s = bsxfun(@minus, data, m);
	v = sum(s.^2, 2);
end %end of subfunction
