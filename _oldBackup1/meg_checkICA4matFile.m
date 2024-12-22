function v = meg_checkICA4matFile (ein, aus, para)
% % % updated 11/07/17 by wp: modify the program for better interface
% % % updated 02/12/16 by wp: modify the program to include MS data
% % % written 20/11/16 by wp: this function is to semi-automatic ICA component rejection

	%% 1. preparison
	% % % some paras
	if nargin < 3 || isempty(para)
		icaWaveDispRange = [-1200 3600];
		skip_id = [];
		freqReductFlag = 1;
		zReductScale = 1;
	else
		icaWaveDispRange = para.icaWaveDispRange;
		skip_id = para.skip_id;
		freqReductFlag = para.freqReductFlag;
		zReductScale = para.zReductScale;
	end
	ListenChar(2);
	% % % indication of starting time
	tmp = clock;
	fprintf('\n===========\nThe ICA semi-manual check start at %04d-%02d-%02d %02d:%02d:%02d.\n', round(tmp));
	% % % get the ica data file structure, but not read in
	try
		fn = matfile(ein);
	catch
		fprintf('error in loading data!');
		return;
	end
	% % % check the file existence
	if exist(aus, 'file')
		fprintf('The previous checking exists, trying to load ...\n');
		try
			tmp = load(aus);
			v = tmp.v;
			v_old = v;
			vFlag = false;
			fprintf(' ...now loaded!\n');
		catch
			fprintf(' ...loading failed!\n');
			vFlag = true;
		end
	else
		vFlag = true;
	end
	% % % construct the file if not there
	if vFlag
		tmp = size(fn, 'ica');
		v.nFs = max(tmp);
		clear tmp;
		v.oc = cell(v.nFs, 1);
		v.label = {'VEOG', 'HEOG', 'EKG', 'LINE', 'MUSCLE', 'JUMP', 'OTHER'};
		v.cf = 1;
		v.cc = 1;
		v_old = [];
	end
	v.ein = ein;
	v.aus = aus;

	%% 2. plot and check
	for is = v.cf : v.nFs
		% % % whether move on
		while KbCheck(-1);
		end
		dataTimer = GetSecs;
		cFlag = 0;
		fprintf('We are going to work on data %02d of %02d, press <space> to continue and <q> to quit.\n', is, v.nFs);
		while ~cFlag
			[tmp1, tmp2, tmp3] = KbCheck(-1);
			if tmp1
				tmp4 = find(tmp3);
				myKey = KbName(tmp4(1));
				switch myKey
					case 'q'
						v.cf = is;
						v.cc = 1;
						save(aus, 'v', 'v_old');
						ListenChar;
						return;
					case 'space'
						cFlag = 1;
					otherwise
						fprintf('\n');
				end
			end
		end
		% % % load ica data
		tmp = fn.ica(is, :);
		if numel(tmp) > 1
			fprintf('The data size was not as expected, please check!');
			return;
		elseif freqReductFlag
			ic = tmp{1};
			x0 = 0.5 : 0.1 : 100;
			[tmp1, tmp2] = min(abs(bsxfun(@minus, ic.xfs', x0)));
			x = ic.xfs(tmp2);
			tmp3 = tmp2(2 : end);
			tmp4 = [tmp2(1 : end - 1); tmp3];
			tmp5 = round(mean(tmp4));
			y = ic.yfs(tmp2, :);
			for k = 2 : length(x) - 1
				y(k, :) = mean(ic.yfs(tmp5(k - 1) : tmp5(k), :));
			end
			ic.xfs = x;
			ic.yfs = y;
			clear tmp1 tmp2 tmp3 tmp5 tmp4 x y x0;
			fprintf('ICA data loaded and power data reduced!\n');
		else
			ic = tmp{1};
			fprintf('ICA data loaded and no reduction was made!\n');
		end

		% % % initialize the output data if not yet
		oc = v.oc{is};
		if ~isfield(oc, 'flag')
			oc.n = ic.nICs;
			if isfield(ic, 'fname')
				oc.fname = ic.fname;
			elseif isfield(ic, 'iPath')
				oc.fname = ic.iPath;
			else
				error('unknown file ID!');
			end
			oc.infoRatio = ic.infoRatio;
			oc.flag = zeros(oc.n, 1, 'single'); %default space
		else
			fprintf('The data was checked partially before, trying to resume...\n');
			if isfield(ic, 'fname') && strcmp(oc.fname, ic.fname)
				fprintf('consistency check succeed!\n');
			elseif isfield(ic, 'iPath') && strcmp(oc.fname, ic.iPath)
				fprintf('consistency check succeed!\n');
			else
				fprintf('consistency check failed!\n');
				fprintf('previous input file name %s, \n', oc.fname);
				fprintf('loaded input file name %s, \n', ic.fname);
				return;
			end
		end
		% % % plot the component templates
		myIn = {'fColor', [1 1 1], 'cvRangeX', icaWaveDispRange, ...
			'cvTmpPtsNum', ic.nxts, 'cvNum', size(ic.yts, 2), 'tpRangeXY', 0.5};
		fp = icaFigParaSet(myIn{:});
		fp.tpTmpPtsNum = round(fp.tpTmpPtsNum / zReductScale);
		ff = icaFigureTemplate(fp);		% plot template
		set(ff.fid, 'visible', 'on');
		theClock = clock;
		fprintf('ICA template prepared @%02d:%02d:%02.0f!\n\n', theClock(4:6));
		theTimer = GetSecs;

		% % % plot for each ic and check
		if is == v.cf
			ii = v.cc;
		else
			ii = 1;
		end

		while ii <= oc.n
			tic;
			if ismember(oc.flag(ii), skip_id)
				if oc.flag(ii)
					fprintf('According to pre-defined criteria, we would skip IC #%02d as<%s>.\n', ii, v.label{oc.flag(ii)});
				else
					fprintf('According to pre-defined criteria, we would skip IC #%02d as<CLEAN>.\n', ii);
				end
				ii = ii + 1;
				continue;
			end
			% % % plot the ICA component
			plotICA_matfile (ff, fp, ic, ii);
			
			% % % do check
			while KbCheck(-1);
			end
			cFlag = 0;
			tFlag = 1;
			while ~cFlag
				[tmp1, tmp2, tmp3] = KbCheck(-1);
				if tmp1
					tmp4 = find(tmp3);
					myKey = KbName(tmp4(1));
					switch myKey
						case 'q'
							v.cf = is;
							v.cc = ii;
							v.oc{is} = oc;
							save(aus, 'v', 'v_old');
							ListenChar;
							return;
						case 'b'
							ii = ii - 2;
							cFlag = 1;
							tFlag = 0;
							if ii + 1 > 0
								fprintf('We would go back to IC #%02d!\n', ii + 1);
							else
								ii = 0;
								fprintf('We could not go back to stone age, no go to IC #%02d!\n', ii + 1);
							end
						case 'v'
							oc.flag(ii) = 1;
							cFlag = 1;
						case 'h'
							oc.flag(ii) = 2;
							cFlag = 1;
						case 'k'
							oc.flag(ii) = 3;
							cFlag = 1;
						case 'l'
							oc.flag(ii) = 4;
							cFlag = 1;
						case 'm'
							oc.flag(ii) = 5;
							cFlag = 1;
						case 'j'
							oc.flag(ii) = 6;
							cFlag = 1;
						case 'n'
							oc.flag(ii) = 7;
							cFlag = 1;
						case 'space'
							oc.flag(ii) = 0;
							cFlag = 1;
						otherwise
							fprintf('\n');
					end
				end
			end
			if tFlag && oc.flag(ii)
				fprintf('IC #%02d of %d was labeled <%s> with %6.2f seconds!\n', ii, ic.nICs, v.label{oc.flag(ii)}, GetSecs - theTimer);
			elseif tFlag
				fprintf('IC #%02d of %d was labeled <CLEAN> with %6.2f seconds!\n', ii, ic.nICs, GetSecs - theTimer);
			end
			ii = ii + 1;
			theTimer = GetSecs;
		end % end of all components
		fprintf('The data %s was processed with %6.2f minutes!\n', oc.fname(end-22:end-8), (GetSecs - dataTimer) / 60);
		v.oc{is} = oc;
		close all;
		save(aus, 'v', 'v_old');
	end % end of all data files

	%% 3. save the data
	v.cf = is;
	v.cc = ii - 1;
	save(aus, 'v', 'v_old');
	ListenChar;
end %end of function


%% sub function 1
function fp = icaFigParaSet(varargin)

	%% 1. para for figure and axis positions etc
	% % % figure properties	
	fp.fColor = [1 1 1];
	fp.fOri = 'Portrait';
% 	fp.fOri = 'Landscape';
% % % 	% % % figure size, about A4 ratio in 1080p screen	
% % % 	fp.fStartXPx = 50;
% % % 	fp.fStartYPx = 50;
% % % 	fp.fWidthPx = 1200;
% % % 	fp.fHeightPx = 850;

	% % % 	figure size, keep A4 ratio for printing
	fp.fgSz = round([0 0 11.7 8.27] * 300);
	
	% % % border settings	
	fp.fStartX = 0.03;	%whole figure start point in X
	fp.fEndX = 0.97;
	fp.fStartY = 0.03;
	fp.fEndY = 0.97;
	fp.bMidX1 = 0.35;	%border of left and right part in X
	fp.bMidX2 = 0.55;	%border of left and right part in X
	fp.bTopY = 0.9;
	fp.bBotY = 0.5;
	fp.spAX1D = 0.03;
	
	%% 2. para for text
	fp.txColorAX	= [1 1 1];
	fp.txRangeX = [0, 1];
	fp.txRangeY = [0, 1];
	fp.txStartX = 0;
	fp.txStartY = 0.5;
	fp.txS = '\n[Here is the Position for component name tag]';
% 	fp.txS2 = '[Here is the Position for correlations]';	
	fp.txSZ = 12;
	fp.txFT = 'Helvetica';
	
	%% 3. curve/trl plot para
	fp.cvColorAX	= [1 1 1];
	fp.cvNum = 10;
	fp.cvRangeX = [0, 1];
	fp.cvRangeY = [0, 10];
	fp.cvTmpPtsNum = 2401;
	fp.cvColor = [0 0 1];
	fp.cvWidth = 0.1;
	fp.trColorAX = [1 1 1];
	
	%% 4. topo plot
	fp.tpColorAX	= [1 1 1];
	fp.tpTmpPtsNum = 1000;
	fp.tpRangeXY = 0.5;
	fp.tpAspRatio = [255 255 1];
	
	%% 5. power plot
	fp.pwColorAX	= [1 1 1];
	fp.pwTmpPtsNum = 100;
	fp.pwRangeX = [0, 200];
	fp.pwRangeY = [0, 1];
	fp.pwColor = [1 0 0];
% 	fp.pwWidth = 2;
	
	%% 6. use iput para to replace default para
	if (rem(length(varargin), 2))
		error('Optional parameters should always go by pairs');
	else
		fs = fieldnames(fp);
		for k = 1 : 2 : (length(varargin)-1)
			if ~ischar (varargin{k}),
				error (['Unknown type of optional parameter name (parameter' ...
				' names must be strings).']);
			else
				yesFieldID = strcmpi(varargin{k}, fs);
				if any(yesFieldID)
					cmd = ['fp.', fs{yesFieldID}, '= varargin{k + 1};'];
					eval(cmd);
				else
					error ('Unknown field names, check again!');
				end
			end
		% change the value of parameter
		end
	end
	
	%% 7. induced parameters
	% % % axes positions
	fp.spAXPos = [fp.spAX1D, fp.spAX1D, -fp.spAX1D, -fp.spAX1D];
	fp.txPos	= [fp.fStartX, fp.bTopY, fp.fEndX - fp.fStartX, fp.fEndY - fp.bTopY] + fp.spAXPos; %text axes position
	fp.cvPos	= [fp.fStartX, fp.fStartY, fp.fEndX - fp.fStartX, fp.bBotY - fp.fStartY] + fp.spAXPos; %curve plot axes position
	fp.tpPos	= [fp.fStartX, fp.bBotY, fp.bMidX1 - fp.fStartX, fp.bTopY - fp.bBotY] + fp.spAXPos; %topo plot axes position
	fp.pwPos	= [fp.bMidX2, fp.bBotY, fp.fEndX - fp.bMidX2, fp.bTopY - fp.bBotY] + fp.spAXPos; %power plot axes position
	fp.trPos	= [fp.bMidX1, fp.bBotY, fp.bMidX2 - fp.bMidX1, fp.bTopY - fp.bBotY] + fp.spAXPos; %trl plot axes position
	
	% % % curve data induction
	fp.cvDataX = linspace(fp.cvRangeX(1), fp.cvRangeX(2), fp.cvTmpPtsNum)';
	tmp = linspace(fp.cvRangeY(1), fp.cvRangeY(2), fp.cvNum + 1);
	tmp(end) = [];
	fp.cvDataY = bsxfun(@plus, (1 + sin(bsxfun(@plus, fp.cvDataX, rand(1, fp.cvNum)) * 20 * pi)) / (2 * fp.cvNum), tmp);
	fp.cvXTick = fp.cvDataX([1, end]);
	fp.cvXTickLabel = num2str(fp.cvXTick);
	fp.cvYTick = linspace(fp.cvRangeY(1), fp.cvRangeY(end), fp.cvNum + 1);
% 	fp.cvYTickLabel = num2str(fp.cvYTick');
	fp.cvYTickLabel = [];

	% % % topo data induction
	fp.tpXY = linspace(-fp.tpRangeXY, fp.tpRangeXY, fp.tpTmpPtsNum);
	fp.tpRangeX = [-fp.tpRangeXY, fp.tpRangeXY];
	fp.tpRangeY = [-fp.tpRangeXY, fp.tpRangeXY];
	% % % 	[fp.tpX, fp.tpY] = meshgrid(tmp, tmp);
	% % % 	fp.tpZ = 1 ./ (sqrt((fp.tpX - fp.tpTmpPtsNum/2).^2 + (fp.tpY - fp.tpTmpPtsNum/2).^2)/fp.tpTmpPtsNum + 1) - 0.5;
% 	fp.tpZ = (fp.tpX - 0.5).^2 + (fp.tpY - 0.5).^2;
	fp.tpZ = bsxfun(@plus, fp.tpXY.^2, (fp.tpXY').^2);

	% % % power data induction
	fp.pwX = linspace(fp.pwRangeX(1), fp.pwRangeX(2), 100)';
	fp.pwY = rand(fp.pwTmpPtsNum, 1) / 3 + (sin(fp.pwX * 4 * pi) + 1) / 3;
	
	% % % trl data
	NN = 1000;
	fp.trx = rand(NN, 1);
	fp.try = 1 : NN;
	
end % end of function

%% sub function 2
function f = icaFigureTemplate(fp)
% % % 28/07/2014 written by wp
% % % plot template for ICA results

	%% 1. figures and axes generation
	f.fid = figure('position', fp.fgSz, 'Units', 'pixels', ...
		'color', fp.fColor, 'visible','off', 'PaperType', 'A4','Resize','off');
	
% % % 	f.fid = figure('position', [fp.fStartXPx fp.fStartYPx fp.fWidthPx fp.fHeightPx], 'Units', 'pixels', ...
% % % 		'color', fp.fColor, 'visible','off', 'PaperType', 'A4','Resize','off');
% % % 	f.fid = figure('position', [0, 0, 8.2677, 11.6929], 'Units', 'inches', ...
% % % 		'color', fp.fColor, 'visible','off', 'PaperType', 'A4','Resize','off');
% % % 	f.fid = figure('Units', 'Normalized', 'color', fp.fColor, ...
% % % 		'visible','off', 'PaperType', 'A4','Resize','off');
% % % 	f.fid = figure('position', [fp.fStartXPx fp.fStartYPx fp.fWidthPx fp.fHeightPx], 'Units', 'Pixels', ...
% % % 		'color', fp.fColor, 'visible','off', 'PaperType', 'A4','Resize','off');
% % % 	'PaperOrientation', fp.fOri, 
	f.txAX = axes('position', fp.txPos, 'Color', fp.txColorAX, 'parent', f.fid);
	f.cvAX = axes('position', fp.cvPos, 'Color', fp.cvColorAX, 'parent', f.fid);
	f.tpAX = axes('position', fp.tpPos, 'Color', fp.tpColorAX, 'parent', f.fid);
	f.pwAX = axes('position', fp.pwPos, 'Color', fp.pwColorAX, 'parent', f.fid);
	f.trAX = axes('position', fp.trPos, 'Color', fp.trColorAX, 'parent', f.fid);
	
	%% 2. text plot
	set(f.txAX, 'xlim', fp.txRangeX, 'ylim', fp.txRangeY);
	f.tx = text(fp.txStartX, fp.txStartY, fp.txS, 'parent', f.txAX, ...
		'fontsize', fp.txSZ, 'fontname', fp.txFT, 'fontweight', 'bold', 'Interpreter', 'none' );
	% % % 	f.tx2 = text(fp.txStartX, fp.txStartY2, fp.txS2, 'parent', f.txAX, ...
	% % % 		'fontsize', fp.txSZ, 'fontname', fp.txFT, 'fontweight', 'bold', 'Interpreter', 'none' );
	set(f.txAX, 'box', 'off', 'XTick', [],'YTick', [], 'visible', 'off');
	
	%% 3. curve line plot
	set(f.cvAX, 'xlim', fp.cvRangeX, 'ylim', fp.cvRangeY);
	f.cvs = plot(fp.cvDataX, fp.cvDataY, 'color', fp.cvColor, 'linewidth', fp.cvWidth, 'parent', f.cvAX);
	set(f.cvAX, 'XTick', fp.cvXTick, 'XTickLabel', fp.cvXTickLabel, 'YTick', fp.cvYTick, 'YTickLabel', fp.cvYTickLabel);
	
	%% 4. topo plot
	set(f.tpAX, 'xlim', fp.tpRangeX, 'ylim', fp.tpRangeY);
% 	f.tp = surface(fp.tpX, fp.tpY, zeros(size(fp.tpZ)), fp.tpZ, 'edgecolor','none', 'parent', f.tpAX);
% 	shading(f.tpAX, 'interp');
	pZ = plt_scale2color(fp.tpZ, 0);
	f.tp = image(fp.tpXY, fp.tpXY, pZ, 'parent', f.tpAX);
	set(f.tpAX, 'DataAspectRatio', fp.tpAspRatio, 'YDir', 'normal', ...
		'XTick', [],'YTick', [], 'box', 'off', 'visible', 'off');
% 	set(f.tpAX, 'box', 'off', 'XTick', [],'YTick', [], 'visible', 'off');
	
	%% 5. power plot
	set(f.pwAX, 'xlim', fp.pwRangeX, 'ylim', fp.pwRangeY);
	f.pw = bar(f.pwAX, fp.pwX, fp.pwY);
	set(f.pw, 'EdgeColor', fp.pwColor, 'FaceColor', fp.pwColor);
% 	f.pw = plot(fp.pwX, fp.pwY, 'color', fp.pwColor, 'linewidth', fp.pwWidth);
	set(f.pwAX, 'box', 'off');
% 	set(f.pwAX, 'box', 'off', 'XTick', [],'YTick', [], 'visible', 'off');

	%% 6. trl plot
	f.tr = plot(fp.trx, fp.try, 'r.', 'parent', f.trAX);
	
end %end of function 


%% sub function 3
function plotICA_matfile (ff, fp, ic, ii)
% % % updated 11/07/17 by wp: this function is now integrated in
% % % written 03/02/17 by wp: this function is to ease plt21_checkICA_matfile

	% % % the text for correlation
	str = [ic.iTag, '_IC', num2str(ii, '%03d')];
	str = sprintf(['<',str, '> r_VEOG = %4.2f, r_HEOG = %4.2f, r_REOG = %4.2f, r_EKG = %4.2f\t'], ic.acc(ii, [1 2 4 3]));
	set(ff.tx, 'String', str);
	
	% % % curves plot
	for k = 1 : size(ic.yts, 2)
		set(ff.cvs(k), 'YData', ic.yts(:, k, ii)/2 + k);
		if abs(ic.xts(1)) < 10 && abs(ic.xts(end)) < 10	%in seconds
			set(ff.cvs(k), 'XData', ic.xts * 1000);
		else
			set(ff.cvs(k), 'XData', ic.xts);	%in ms
		end					
	end
	tmpX = [-1 0 1 2 3] * 1000;
	tmpY = 1 : size(ic.yts, 2);
	set(ff.cvAX, 'YLim', [0 size(ic.yts, 2) + 1], 'fontsize', 20, ...
	'XTIck', tmpX', 'XtickLabel', num2str(tmpX'), 'xlim', tmpX([1, end]), ...
	'YTIck', tmpY, 'YtickLabel', num2str(tmpY'));

	% % % topo plot
	pZ = plt_topo2Data(ic.Z(:, ii), ic.loc2d, fp.tpRangeXY, fp.tpTmpPtsNum);
	sid = isnan(pZ);
	pZ = plt_scale2color(pZ, 0);
	for k = 1 : 3
		tmp = pZ(:, :, k);
		tmp(sid) = 1;
		pZ(:, :, k) = tmp;
	end
	set(ff.tp, 'CData', pZ);

	% % % power plot		
	set(ff.pw, 'XData', ic.xfs, 'YData', ic.yfs(:, ii));
	set(ff.pwAX, 'XLim', [0 ic.fM], 'XTick', [2 4 10 30 50 ic.fM], 'XTickLabel', num2str([2 4 10 30 50 ic.fM]'));

	% % % trl plot
	set(ff.tr, 'XData', ic.trv(:, ii), 'YData', 1 : ic.nm);
	set(ff.trAX, 'XLim', [0 1], 'XTIck', [0 0.5 1], 'YLim', [0 ic.nm + 1]);
	drawnow;
	
end %end of function