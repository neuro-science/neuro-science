	function [t, p]= mat_tThreshFromP(n, p0, delta, tStart, twoSideFlag)
	% % % 22/06/17 written by wp

		%% 1. check paras
		if nargin < 1
			fprintf('You need input the number of cases!\n')
			return;
		end
		
		if nargin < 2 || isempty(p0)
			p0 = 0.05;
		end
		
		if nargin < 3 || isempty(delta)
			delta = 1e-5;
		end
		
		if nargin < 4 || isempty(tStart)
			tStart = 2.2;
		end
		
		if nargin < 5 || isempty(twoSideFlag)
			twoSideFlag = true;
		end
		
		scale = 0.1;
		if twoSideFlag
			x = 2;
		else
			x = 1;
		end
		
		%% 2. do the recursive loop
		goFlag = 1;
		t = tStart;
		while goFlag
			p = x * (1 - tcdf(t, n - 1));
			if p - p0 > delta
				t = t * (1 + scale);
			elseif p - p0 < -delta
				t = t * (1 - scale);
			else
				goFlag = 0;
			end
		end
	end
