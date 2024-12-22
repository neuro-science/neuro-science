function p = plt19_srconsurface2_debug(vc, fc, cv, si, vn, d, varargin)
% % % updated 05/07/2017 by wp: add fid
% % % written 06/01/16 by wp: for dual layer display.
% % % Instead of modifying the anatomy data, we superimpose another layer
% % % on it.


	%% prepare
	% % % default parameters for viewing
 	p.figPosSize = [900 100 900 900]; %position and size
	p.figColor = [1 1 1]; %background of figure
	p.axPosSize = [0 0 1 1];
	p.axLimSetFlag = 1; %axis limit need to be set
	p.axLim = [-100 100];
	p.axColor = [1 1 1]; %background of axis
	p.markerSize = 50;	%size of sources
	p.srcColor = [1 0 0];
	p.baseColorScale = 0.5;
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
	p.barPossible = 0;
	p.alpha4function = 0.75;
	p.functionExtendDistance = 0.2; %2mm above anatomy
	p.curvThresh = 0; %Threshold for curvatures (gyrus/sulcus)
	
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
			case 'sc'
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
			case 'fd'
				p.functionExtendDistance = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'av'
				p.alpha4function = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'vh'
				p.curvThresh = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'fid'
				p.fid = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			otherwise
				nUnknownPara = nUnknownPara + 1;
				fprintf('\nUnknown parameter #%d: " %s " ...\n', varargin{kCounter});
		end
		kCounter = kCounter + 1;
	end % readin end
	
	% % % color settings for anatomy
 	cortex_light = [0.781 0.762 0.664];
	cortex_dark  = [0.781 0.762 0.664]/2;
	if ~isempty(cv)
		cv(cv >= p.curvThresh) = 1;
		cv(cv < p.curvThresh) = -1;
		p.crxColorBack = ((1 + cv) * cortex_dark + (1-cv) * cortex_light)/2;
	else
		p.crxColorBack = repmat(cortex_light, [size(vc, 1), 1]);
	end
	
	% % % color settings for function, with subsample
	if ~isfield(p, 'subsample') || isempty(p.subsample)
		p.subsample = 0.1;
	elseif p.subsample > 1
		p.subsample = 1./p.subsample;
	end
	vc1 = vc + p.functionExtendDistance * vn; %tri1 = fc;
	p.crxColorFront = ones(size(vc1));
	p.vertexAlpha = zeros(size(vc1, 1), 1);

	% % % para set from fixed input
	p.srcLocations = vc1(si, :);
	if isempty(d)
		fprintf('Only location provided, uniform color [default: red] will be used.\n');
		p.srcPlot = 1;
	else
		p.barPossible = 1;
		fprintf('Both location and strength provided, color illustrate strength.\n');
		p.srcPlot = 0;
		switch p.intMethod
			case 'gaussian'
				vout = plt07_gaussian_dist(vc1, p.srcLocations, p.gaussFilterThresh, d);
			case 'nearest'
				vout = interp_ungridded(p.srcLocations, vc1, 'data', d, 'projmethod', 'nearest');
			case 'sphere_avg'
				vout = interp_ungridded(p.srcLocations, vc1, 'data', d, 'projmethod', 'sphere_avg', 'sphereradius', p.gaussFilterThresh);
			case 'sphere_weighteddistance'
				vout = interp_ungridded(p.srcLocations, vc1, 'data', d, 'projmethod', 'sphere_weighteddistance', 'sphereradius', p.gaussFilterThresh);
			case 'smudge'
				vout = interp_ungridded(p.srcLocations, vc1, 'data', d, 'projmethod', 'smudge');
			otherwise
				error('unknown projection method!\n');
		end
		[cTmp, sid, theMap, cMax] = plt03_scale2color(vout, p.colorThresh, p.colorScale, p.colorStyle);
		p.crxColorFront(sid, :) = cTmp(sid, :);
		p.vertexAlpha(sid) = p.alpha4function;
		txid = find(sid(si, :) > 0);
		clear cTmp sid;
	end
	
	%% show images
	% % % prepare figue
	% % % figure size etc
	if ~isfield(p, 'fid')
		p.fid = figure;
	end
	% % % render
	set(p.fid, 'position', p.figPosSize, 'Color', p.figColor, 'renderer', 'OpenGL');
	% % % 	loop for 3 axes
	p.ax = axes('position', p.axPosSize, 'parent', p.fid, 'Color', p.axColor, ...
	 'XTick', [], 'XLim', p.axLim, 'YTick', [], 'YLim', p.axLim, 'ZTick', [], 'ZLim', p.axLim); %create axis
	hold on;
	p.mH = patch('faces', fc, 'vertices', vc, 'parent', p.ax, ...
		'facecolor', 'interp',  'edgecolor', 'none', ...
		'facealpha', 0.99, 'FaceVertexCData', p.crxColorBack);
	p.mH1 = patch('faces', fc, 'vertices', vc1, 'parent', p.ax, ...
		'facecolor', 'interp',  'edgecolor', 'none', ...
		'FaceVertexCData', p.crxColorFront, 'AlphaDataMapping', 'none',...
		'FaceVertexAlphaData', p.vertexAlpha, 'facealpha', 'interp');
	if p.srcPlot
		p.dotLocations = p.srcLocations + bsxfun(@times, abs(p.srcLocations), p.dotExtRatio * p.viewDir);
% 		p.dotLocations = bsxfun(@times,	p.srcLocations, (1 + p.dotExtRatio));
		p.dH = plot3(p.dotLocations(:, 1), p.dotLocations(:, 2), p.dotLocations(:, 3), '.', ...
			'markersize', p.markerSize, 'color', p.srcColor, 'parent', p.ax);
	end
	if isfield(p, 'txtSubscript')
		p.txtLocations = p.srcLocations + bsxfun(@times, abs(p.srcLocations), p.txtExtRatio * p.viewDir);
		if isnumeric(p.txtSubscript)
			tmp = num2str(p.txtSubscript');
	% 			p.txtLocations = bsxfun(@minus, p.srcLocations, (1 + p.txtExtRatio));
			p.txtLocations = p.txtLocations(p.txtSubscript, :);
			p.dT = text(p.txtLocations(:, 1), p.txtLocations(:, 2), p.txtLocations(:, 3), tmp, 'fontsize', 20);
			clear tmp;
		elseif strcmp(p.txtSubscript, 'act')
			p.txtLocations = p.txtLocations(txid, :); %bsxfun(@minus, p.txtLocations(txid, :), [1 -5 0]);
			for ttt = 1 : length(txid)
% 				p.dT(ttt) = text(p.txtLocations(ttt, 1), p.txtLocations(ttt, 2), p.txtLocations(ttt, 3), ...
% 					[num2str(d(txid(ttt)), '%0.1f'), '@', num2str(txid(ttt))], ...
% 					'fontsize', 20, 'HorizontalAlignment', 'center', 'Fontweight', 'bold');
				p.dT(ttt) = text(p.txtLocations(ttt, 1), p.txtLocations(ttt, 2), p.txtLocations(ttt, 3), ...
					num2str(txid(ttt)), ...
					'fontsize', 10, 'HorizontalAlignment', 'center', 'Fontweight', 'bold');
			end
			clear tmp;
		elseif strcmp(p.txtSubscript, 'loc') %%%buggy
			for ttt = 1 : length(d)
				p.dT(ttt) = text(p.txtLocations(ttt, 1), p.txtLocations(ttt, 2), p.txtLocations(ttt, 3), ...
					{[num2str(ttt), '__', num2str(ttt)]}, ...
					'fontsize', 20, 'color', [0 1 0], 'HorizontalAlignment', 'center', 'Fontweight', 'bold');
			end
			clear tmp;
		else
			fprintf('\nUnknown text parameter #%d: " %s " ...\n', p.txtSubscript);
		end
	end
	if ~isempty(p.bar) && p.barPossible
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
	p.vc1 = vc1;p.vc = vc;
end

