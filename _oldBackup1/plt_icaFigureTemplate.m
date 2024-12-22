function f = plt_icaFigureTemplate(fp)
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

