function c = plt01_genColorRGB( N0, M, pltFlag)
% % % 06/07/14	written by wp to generate RGB colors
% % % N0 more colors will be generated.  
% % % default is 0-1 in each channel, or integer scaled to M

	%% check inputs
	if nargin < 3
		pltFlag = 0;
	end
	if nargin < 2 || M == 1
		intFlag = 0;
	else
		intFlag = 1;
	end
	N = max(ceil((N0 + 2)^(1/3)), 2);
	
	%% do it
	x = linspace(0, 1, N);
	c = zeros(N^3, 3);
	ct = 0;
	for k1 = 1 : N
		for k2 = 1 : N
			for k3 = 1 : N
				ct = ct + 1;
				c(ct, 1) = x(k1);
				c(ct, 2) = x(k2);
				c(ct, 3) = x(k3);
			end
		end
	end
	c([1, end], :) = [];
	
	%% sort
	[y, I] = sort(std(c, 0, 2), 'descend');
	c = c(I, :);
	
	%% plot if necessary
	if pltFlag
		npts = 100;
		y = linspace(0, 10, npts);
		figure;
		hold on;
		for k = 1 : size(c, 1)
			x = zeros(npts, 1) + k;
			plot(x, y, 'color', c(k, :), 'linewidth', 2);
		end
	end
	
	%% modify for out put
	if intFlag
		c = round(c * M);
	end

	
end

