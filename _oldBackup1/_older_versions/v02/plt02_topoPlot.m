function cfg = plt02_topoPlot(data, loc2d, cfg)
% % % updated 05/08/2014 - make it more general
% % % updated 12/06/2014 - locations as input para rather than fields
% % % updated 16/05/2013 - improve figure property
% % % updated 15/05/2013 - add central sensor bold
% % % created by wp for sensor topo plot

% % % data  :	nchans x 1, data
% % % loc2d :	nchans x 2	locations 2d (x, y)
% % % cfg   :	parameters for plotting


	%% headers
	% % % check input
	if nargin < 2
		error('at least data and locations shall be provided!');
	elseif nargin < 3
		cfg = [];
	end
	
	% % % set some parameters
	if ~isfield(cfg, 'exRange')
		cfg.exRange = 0.5;
	end
	if ~isfield(cfg, 'nPixels1D')
			cfg.nPixels1D = 1000;
	end
	if ~isfield(cfg, 'nPixels1D')
			cfg.nPixels1D = 1000;
	end
	if ~isfield(cfg, 'doPlot')
			cfg.doPlot = 1;
	end
	if ~isfield(cfg, 'max')
			cfg.max = [];
	end

	% % % prepare data
	if verLessThan('matlab', '8.1')
		F = TriScatteredInterp(loc2d, data);
	else
		F = scatteredInterpolant(loc2d, data);
	end
	tmp = linspace(-cfg.exRange, cfg.exRange, cfg.nPixels1D);
	[pX,pY] = meshgrid(tmp, tmp);
	pZ = F(pX, pY);

	%% Plot
	if cfg.doPlot
		% % % figure handles
		if ~isfield(cfg, 'fid')
			cfg.fid = figure('position', [900 100 900 900], 'Color', [1 1 1]);
			cfg.aid = axes('position', [0.1 0.1 0.8 0.8]);
		elseif ~isfield(cfg, 'aid')
			cfg.aid = axes('position', [0.1 0.1 0.8 0.8]);
		end
		
		% % % color 
		cZ = plt03_scale2color(pZ(:), 0, cfg.max);
		
		% % % Remove outsiders
		cZ = reshape(cZ, [cfg.nPixels1D * cfg.nPixels1D, 3]);
		outSiders = (sqrt(pX.^2 + pY.^2) > cfg.exRange);
		nanZ = isnan(cZ(:, 1)) | outSiders(:);
		cZ(nanZ, :) = 1;
		cZ = reshape(cZ, [cfg.nPixels1D, cfg.nPixels1D, 3]);
		
		% % % plot		
		h = image(pX(1, :), pY(:, 1), cZ, 'parent', cfg.aid);
		set(gca, 'ydir', 'normal');
		
		% % % landmark electrodes on?
		if ~isfield(cfg, 'elecOn')
				cfg.elecOn = 0;
		elseif cfg.elecOn 
			if ~isfield(cfg, 'elecSize')
				cfg.elecSize = 10;
			end
			if isfield(cfg, 'spSensors')
				cfg.spFlag = 1;
				if ~isfield(cfg, 'elecSizeB')
					cfg.elecSizeB = 20;
				end
			end
		end
		
		% % % plot electrodes if needed
		if cfg.elecOn
			hold on;
			h2 = plot(loc2d(:, 1), loc2d(:, 2), '.k', 'markersize', cfg.elecSize);
			if isfield(cfg, 'spSensors')
				h3 = plot(loc2d(cfg.spSensors, 1), loc2d(cfg.spSensors, 2), '.m', 'markersize', cfg.elecSizeB);
			end
			if isfield(cfg, 'labels')
				text(loc2d(:, 1), loc2d(:, 2), cfg.labels);
			end
		end
		axis(cfg.aid, [-cfg.exRange cfg.exRange -cfg.exRange cfg.exRange]);
		axis equal;
		axis off;
	else
% 		h = get(cfg.aid, 'children');
% 		set(h(end), 'CData', pZ);
	end
	
end %end of function