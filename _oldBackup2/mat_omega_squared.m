function w2 = mat_omega_squared(X, ids)

% % % written on 26/03/2018 by wp, compute omega squared
% % % inputs:
% % %		X (n x m), n is cases
% % %		ids (n x 1), indices for conditions
% % % outputs:
% % %		w2 - omega sqaured

	% % % check inputs
	[n, m] = size(X);
	if size(ids, 1) ~= n
		error('data size mismatch!');
	end
	
	% % % sort indices
	[grp, idx] = sort(ids);
	X = X(idx, :);
	[x1, x2] = unique(grp);
	N = length(x1);
	idN = [x2; n+1];
	nS = idN(2 : end) - idN(1 : end - 1);
	
	% % % compute
	meanGrand = mean(X, 1);
	meanGroup = zeros(N, m);
	ssErrTmp = zeros(N, m);
	for k = 1 : N
		meanGroup(k, :) = mean(X(idN(k) : idN(k+1) - 1, :), 1);
		ssErrTmp(k, :) = sum(bsxfun(@minus, X(idN(k) : idN(k+1) - 1, :), meanGroup(k, :)).^2, 1);
	end
	ssGroup = nS' * (bsxfun(@minus, meanGroup, meanGrand).^2);
	ssErr = sum(ssErrTmp, 1);
	% within-groups MS
	msErr = ssErr / (n - N);
	% SS_t, the sum of both
	ssTot = ssGroup + ssErr;
	% corresponding df
	w2 = (ssGroup - (N - 1) * msErr) ./ (ssTot + msErr);
	
end


	%end of function