function cfg = plt_tfPlot2d( data, cfg )
% % % 22/05/24 modified by wp
% % % low level improvement for newer colormap

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

	%% 1. read data information and prepare environment
	% % % get information of data, the last dimension 3 for rgb, 1 for unicolor
	[ny, nx, nc] = size(data);	
	if nc ~= 1
		fprintf('The input data is not 2D: %d pages!\n', nc);
	end
	
	%% 2. fast plot of certain data without input parameters
	if nargin < 2
		% % % convert the data into color via all default settings		
		[theColor, alpha] = plt_scale2color4(data);

		% % % xy range		
		x = 1 : nx;
		y = 1 : ny;
		
		% % % plot
% 		dispColor = bsxfun(@times, theColor, alpha);
% 		cfg.h = image('XData', x, 'YData', y, 'CData', dispColor);
		cfg.h = image('XData', x, 'YData', y, 'CData', theColor);
		set(gca, 'ydir', 'normal');
		return;
	end

	%% 3. prepare those with paras, give defaults when empty
	% % % set figure
	if ~isfield(cfg, 'fid')
		cfg.fid = figure('Position', [500 100 1000 1000], 'PaperType', 'A4', 'Color', [1 1 1]);
	end
	
	% % % set ax
	if ~isfield(cfg, 'aid')
		cfg.aid = axes('position', [0.2 0.2 0.7 0.7], 'parent', cfg.fid);
	end
	
	% % % color threshold	
	if ~isfield(cfg, 'th')
		cfg.th = 0;
	end
	
	% % % colormap
	if ~isfield(cfg, 'colormap')
		cfg.colormap = 'm';
	end
	
	%% 4. do plotting
	% % % 	whether interpolate
	if isfield(cfg, 'szInterpolate')	&& ~isempty(cfg.szInterpolate) %do interpolation
		
		[y, x] = ind2sub([ny, nx], 1 : nx * ny);
		
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
			[cZ, cAlpha] = plt_scale2color4(pZ, cfg.th, cfg.cLim, cfg.colormap);
		else
			[cZ, cAlpha] = plt_scale2color4(pZ, cfg.th, [], cfg.colormap);
		end

% 		cZ = reshape(cZ, [size(pZ), 3]);
% 		cfg.h = image(dX, dY, cZ, 'parent', cfg.aid);
		cfg.h = image('XData', dX, 'YData', dY, 'CData', cZ, 'parent', cfg.aid);
		
	else %plot the original data point

		% % % get the color as triplets
		if isfield(cfg, 'cLim')
			cMap = plt_scale2color4(data, cfg.th, cfg.cLim, cfg.colormap);
		else
			cMap = plt_scale2color4(data, cfg.th, [], cfg.colormap);
		end
		
		% % % do plot
		cfg.h = image('XData', 1 : nx, 'YData', 1 : ny, 'CData', cMap, 'parent', cfg.aid);
	end
	
	%% 5. set axis and others
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

