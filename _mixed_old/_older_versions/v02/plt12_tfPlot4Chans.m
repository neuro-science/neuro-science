function [pZ, pX, pY] = plt12_tfPlot4Chans(data, posiChan, cfg)
% % % written 25/11/2014

	%% 1. para set
	% % % check fig parameters
	if nargin < 3 || ~isfield(cfg, 'fid') || isempty(cfg.fid)
		cfg.fid = figure('position', [900 100 1456 900], 'Color', [1 1 1]);
	end
	
	% % % check chan parameters
	if nargin < 3 || ~isfield(cfg, 'ratio') || isempty(cfg.ratio)
		cfg.ratio = 0.8;
	end
	
	% % % check chan parameters
	if nargin < 3 || ~isfield(cfg, 'size') || isempty(cfg.size)
		cfg.size = [0.1, 0.1];
	end
	
	% % % check chan parameters
	if nargin < 3 || ~isfield(cfg, 'a_off') || isempty(cfg.a_off)
		cfg.a_off = 1;
	end
	
	% % % check nch
	nCh1 = size(data, 3);
	nCh2 = size(posiChan, 1);
	if nCh1 == nCh2
		n = nCh1;
	else
		error('channel numbers mismatch!');
	end
	pos = normalLocs(posiChan, cfg.ratio);
	
	%% 2. plot
	for k = 1 : n
		cfg.aid = axes('position', [pos(k, 1), pos(k, 2), cfg.size(1), cfg.size(2)], 'parent', cfg.fid);
		plt08_tfPlot(data(:, :, k), cfg);
		if cfg.a_off
			set(cfg.aid, 'visible', 'off');
		end
	end
end %end of function


function y = normalLocs(x, r)
	y = zeros(size(x));
	for k = 1 : size(x, 2)
		z = x(:, k);
		y(:, k) = r * (z - min(z)) / (max(z) - min(z)) + (1-r)/2; 
	end
end