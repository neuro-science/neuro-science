function cfg = plt_tfPlot2d( data, cfg )
% % % 21/01/19 modified by wp
% % % fast plot slight improvement

% % % 20/12/17 modified by wp? for the new naming rule

% % % 22/08/14 modified by wp
% % % more general

% % % 19/05/14 modified by wp
% % % more margins

% % % 13/05/14 modified by wp
% % % possible empty szInterpolate

% % % 12/05/14 re-written by wp
% % % plot data as 2-D color image

	%% read data information and prepare environment
	% % % get information of data
	[ny, nx] = size(data);
	
	%% fast plot of certain data if no parameters
	if nargin < 2
		% % % set limits
		cfg.maxLimColor = max(abs(data(:)));
		% % % transformation
		data = 64 * data ./ cfg.maxLimColor;
		% % % xy range		
		x = 1 : nx;
		y = 1 : ny;
		% % % plot
		h = image(x, y, data);
		set(gca, 'ydir', 'normal');
		return;
	end

	% % % set figure
	if ~isfield(cfg, 'fid')
		cfg.fid = figure('Position', [500 100 1000 1000], 'PaperType', 'A4', 'Color', [1 1 1]);
	end
	% % % set ax
	if ~isfield(cfg, 'aid')
		cfg.aid = axes('position', [0.2 0.2 0.7 0.7], 'parent', cfg.fid);
	end
	
	%% more complicated
	% % % 	whether interpolate
	if isfield(cfg, 'szInterpolate')	&& ~isempty(cfg.szInterpolate) %do interpolation
		[y, x] = ind2sub([ny, nx], 1 : nx*ny);
	% % % version check needed		
		if verLessThan('matlab', '8.1')
			F = TriScatteredInterp(x', y', data(:));
		else
			F = scatteredInterpolant(x', y', data(:));
		end
		dX = linspace(1, nx, cfg.szInterpolate * (nx - 1) + 1);
		dY = linspace(1, ny, cfg.szInterpolate * (ny - 1) + 1);
		[pX, pY] = meshgrid(dX, dY);
		pZ = F(pX, pY);
		if isfield(cfg, 'cLim')
			cZ = plt_scale2color(pZ(:), 0, cfg.cLim);
		else
			cZ = plt_scale2color(pZ(:), 0);
		end
		cZ = reshape(cZ, [size(pZ), 3]);
		h = image(dX, dY, cZ, 'parent', cfg.aid);
	else %plot the original data point
		% % % get the color as triplets
		if isfield(cfg, 'cLim')
			cMap = plt_scale2color(data(:), 0, cfg.cLim);
		else
			cMap = plt_scale2color(data(:), 0);
		end
		cData = reshape(cMap, [ny, nx, 3]);
		clear data cMap;
		% % % do plot
		h = image(1 : nx, 1 : ny, cData, 'parent', cfg.aid);
	end
	
	% % % set axis
	if ~isfield(cfg, 'noAxLabel')
		if ~isfield(cfg, 'yPos')
			cfg.yPos = 1 : ny;
		end
		if ~isfield(cfg, 'xPos')
			cfg.xPos = 1 : nx;
		end
		if ~isfield(cfg, 'yAxisLabel')
			cfg.yAxisLabel = num2str((1 : ny)');
		end
		if ~isfield(cfg, 'xAxisLabel')
			cfg.xAxisLabel = num2str((1 : nx)');
		end
		set(cfg.aid, 'DataAspectRatio', [255 255 1], ...
			'ydir', 'normal', 'ytick', cfg.yPos, 'YTickLabel', cfg.yAxisLabel, ...
			'xtick', cfg.xPos, 'XTickLabel', cfg.xAxisLabel, 'TickDir','out');
		if isfield(cfg, 'fontsize')
			set(cfg.aid, 'fontsize', cfg.fontsize);
		end
		if isfield(cfg, 'title')
			title(cfg.title);
		end
		if isfield(cfg, 'yLabel')
			ylabel(cfg.yLabel);
		else
			ylabel('Frequency(Hz)');
		end
		if isfield(cfg, 'xLabel')
			h1 = xlabel(cfg.xLabel);
		else
			xlabel('Time(ms)');
		end
		if isfield(cfg, 'xlim')
			set(cfg.aid, 'xlim', cfg.xlim);
		end
		if isfield(cfg, 'ylim')
			set(cfg.aid, 'ylim', cfg.ylim);
		end
	else
		set(cfg.aid, 'visible', 'off');
		set(cfg.aid, 'DataAspectRatio', [255 255 1], 'ydir', 'normal');
	end
	
	
end %end of function

