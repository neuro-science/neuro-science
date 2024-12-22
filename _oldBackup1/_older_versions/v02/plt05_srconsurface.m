function p = plt05_srconsurface(vc, tri, cv, src, ax, varargin)
% % % modified 05/01/16 by wp: for curvatures, remove subsample.
% % % modified 29/09/14 by wp: for smoother figures.
% % % interp_ungridded.m in fieldtrip was used for some interplation of
% % % data, thanks to the authors!!!


% % % modified 08/08/14 by wp: for mroe general use.
% % % modified 12/05/14 by wp: new thresholds etc.
% % % rewrote  09/05/14 by wp: faces updated as in freesurfer_read_surf suggested

	%% prepare
	% % % default parameters for viewing
 	p.figPosSize = [900 100 900 900]; %position and size
	p.figColor = [1 1 1]; %background of figure
	p.axPosSize = [0 0 1 1];
	p.axLimSetFlag = 1; %axis limit need to be set
	p.axLim = [-100 100];
	p.axColor = [1 1 1]; %background of axis
	p.markerSize = 20;	%size of sources
	p.srcColor = [1 0 0];
	p.baseColorScale = 0.5;
	p.faceAlpha = 1;
	p.gaussFilterThresh = 20;
	p.viewDir = [0, -1, .7];
	p.dotExtRatio = 0;
	p.txtExtRatio = 0;
   p.colorThresh = 0.2;
	p.intMethod = 'gaussian';
	p.colorStyle = 'jet';
	p.colorScale = [];
	p.bar = [];
	p.numIter = 1;
	p.maxSrcNum = 10000;	% assumption of maximum sources
	p.barPosible = 0;

	% % % input parameter overwrite defaults
	nVarArgIn = length(varargin);
	kCounter = 1; %arg in counter
	nUnknownPara = 0;	%number of unkown input, display only
	while kCounter < nVarArgIn
		switch lower(varargin{kCounter})
			case 'figpossize'
				p.figPosSize = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'axpossize'
				p.axPosSize = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'axlim'
				p.axLim = varargin{kCounter + 1};
				p.axLimSetFlag = 1;
				kCounter = kCounter + 1;
			case 'axcolor'
				p.axColor = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'markersize'
				p.markerSize = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'facealpha'
				p.faceAlpha = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'viewdir'
				p.viewDir = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'dr'
				p.dotExtRatio = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'tid'
				p.txtSubscript = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'tr'
				p.txtExtRatio = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'gth'
				p.gaussFilterThresh = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'cth'
				p.colorThresh = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'cs'
				p.colorScale = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'cy'
				p.colorStyle = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'si'
				p.srcID = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'ni'
				p.numIter = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'itm'
				p.intMethod = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'bar'
				p.bar = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			otherwise
				nUnknownPara = nUnknownPara + 1;
				fprintf('\nUnknown parameter #%d: " %s " ...\n', varargin{kCounter});
		end
		kCounter = kCounter + 1;
	end % readin end
	
	% % % color settings
 	cortex_light = [0.781 0.762 0.664];
	cortex_dark  = [0.781 0.762 0.664]/2;
	if ~isempty(cv)
		cv(cv >= 0) = 1;
		cv(cv <= 0) = -1;
		p.crxColor = ((1 + cv) * cortex_light + (1-cv) * cortex_dark)/2;
	else
		p.crxColor = repmat(cortex_light, [size(vc, 1), 1]);
	end
	
	% % % para set from fixed input
	switch size(src, 2) 
		case 3
				fprintf('Only location provided, uniform color [default: red] will be used.\n');
				p.srcPlot = 1;
				p.srcLocations = src;
		case 4
			p.barPosible = 1;
			fprintf('Both location and strength provided, color illustrate strength.\n');
			p.srcPlot = 0;
			switch p.intMethod
				case 'gaussian'
					vout = plt07_gaussian_dist(vc, src(:, 1:3), p.gaussFilterThresh, src(:, 4));
				case 'nearest'
					vout = interp_ungridded(src(:, 1:3), vc, 'data', src(:, 4), 'projmethod', 'nearest');
				case 'sphere_avg'
					vout = interp_ungridded(src(:, 1:3), vc, 'data', src(:, 4), 'projmethod', 'sphere_avg', 'sphereradius', p.gaussFilterThresh);
				case 'sphere_weighteddistance'
					vout = interp_ungridded(src(:, 1:3), vc, 'data', src(:, 4), 'projmethod', 'sphere_weighteddistance', 'sphereradius', p.gaussFilterThresh);
				case 'smudge'
					vout = interp_ungridded(src(:, 1:3), vc, 'data', src(:, 4), 'projmethod', 'smudge');
				otherwise
					error('unknown projection method!\n');
			end
			[cTmp, sid, theMap, cMax] = plt03_scale2color(vout, p.colorThresh, p.colorScale, p.colorStyle);
			p.crxColor(sid, :) = cTmp(sid, :);
			clear cTmp sid;
		case 6
			fprintf('Both location and color provided, individual color will be used.\n');
			p.srcPlot = 1;
			p.srcLocations = src(:, 1:3);
			p.srcColor = src(:, 4:6);
		case 7
			fprintf('Both location and color provided, size as well, individual color will be used.\n');
			p.srcPlot = 1;
			p.srcLocations = src(:, 1:3);
			p.srcColor = src(:, 4:6);
			p.markerSize = src(:, 7);
		otherwise
			fprintf('I accept 3, 4, 6 or 7 columns input, now no source will be shown.\n');
			p.srcPlot = 0;
	end
	
	%% show images
	% % % prepare figue
	if isempty(ax)
		% % % figure size etc
		p.fid = figure('position', p.figPosSize, 'Color', p.figColor);
		% % % render
		set(p.fid, 'renderer', 'OpenGL');
		% % % 	loop for 3 axes
		p.ax = axes('position', p.axPosSize, 'parent', p.fid, 'Color', p.axColor, ...
		 'XTick', [], 'XLim', p.axLim, 'YTick', [], 'YLim', p.axLim, 'ZTick', [], 'ZLim', p.axLim); %create axis
	else
		p.ax = ax;
	end
	hold on;
	p.mH = patch('faces', tri(:, [1 3 2]), 'vertices', vc, 'parent', p.ax, ...
		'facecolor', 'interp',  'edgecolor', 'none', ...
		'facealpha', p.faceAlpha, 'FaceVertexCData', p.crxColor);
	if p.srcPlot
		p.dotLocations = p.srcLocations + bsxfun(@times, abs(p.srcLocations), p.dotExtRatio * p.viewDir);
		p.dH = plot3(p.dotLocations(:, 1), p.dotLocations(:, 2), p.dotLocations(:, 3), '.', ...
			'markersize', p.markerSize, 'color', p.srcColor, 'parent', p.ax);
		if isfield(p, 'txtSubscript')
			tmp = num2str(p.txtSubscript');
			p.txtLocations = p.srcLocations + bsxfun(@times, abs(p.srcLocations), p.txtExtRatio * p.viewDir);
			p.txtLocations = p.txtLocations(p.txtSubscript, :);
			p.dT = text(p.txtLocations(:, 1), p.txtLocations(:, 2), p.txtLocations(:, 3), tmp, 'fontsize', 20);
		end
	end
	if ~isempty(p.bar) && p.barPosible
		colormap(theMap);
		p.cb = colorbar(p.bar);
		set(p.cb, 'ytick', [0, 1], 'yticklabel', num2str(cMax));
	end
	axis tight;
	lighting gouraud;
	view(p.viewDir);
	set(p.ax, 'DataAspectRatio', [1 1 1], 'box', 'off', 'visible', 'off');
% 	p.l = camlight('headlight');
% 	set(p.l, 'color', [.2 .2 .2]'); %.1 .1 .1 , 'Visible', 'off

end

