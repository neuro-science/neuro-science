
function d = cmp38_clusterCoherence(c, method, t, num)

% % % This function summarize the t values inside the clusters
% % % written in 17/10/14 by wp

	%% check inputs
	if nargin < 2
		error('You need to specify the cluster and method!\n');
	elseif nargin < 3
		warning('t-value not provided, use number of connections instead!\n');
		tflag = 0;
	else
		tflag = 1;
	end
	
	%% do for 1 of 3 methods
	switch lower(method)
		% % % time-frequency
		case {'tf'} % time-frequency
			% % % get the ranges
			if nargin < 4 || ~isfield(num, 'F') || isempty(num.F)
				num.F = max(c.ed(:, 3));	%check later
			end
			if nargin < 4 || ~isfield(num, 'T') || isempty(num.T)
				num.T = max(c.ed(:, 4));
			end
			d = zeros(num.F, num.T);
			if tflag
				for k = 1 : c.ne
					d(c.ed(k, 3), c.ed(k, 4)) = d(c.ed(k, 3), c.ed(k, 4)) + t(c.id(k));
				end
			else
				for k = 1 : c.ne
					d(c.ed(k, 3), c.ed(k, 4)) = d(c.ed(k, 3), c.ed(k, 4)) + 1;
				end
			end
		case {'sp'} % spatial
			if nargin < 4 || ~isfield(num, 'V') || isempty(num.V)
				num.V = max(c.nd);
			end
			d = zeros(num.V, 1);
			if tflag
				for k = 1 : c.ne
					d(c.ed(k, 1)) = d(c.ed(k, 1)) + t(c.id(k));
					d(c.ed(k, 2)) = d(c.ed(k, 2)) + t(c.id(k));
				end
			else
				for k = 1 : c.ne
					d(c.ed(k, 1)) = d(c.ed(k, 1)) + 1;
					d(c.ed(k, 2)) = d(c.ed(k, 2)) + 1;
				end
			end
		otherwise
			error('unknown methods');
	end

	
end % end of function