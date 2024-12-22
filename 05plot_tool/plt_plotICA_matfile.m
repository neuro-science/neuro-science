function plt_plotICA_matfile (ff, fp, ic, ii, loc2d)
	if nargin < 5 || isempty(loc2d)
		loc2d = ic.loc2d;
	end
% % % written 03/02/17 by wp: this function is to ease plt21_checkICA_matfile

	% % % the text for correlation
	str = [ic.iTag, '_IC', num2str(ii, '%03d')];
	str = sprintf(['<',str, '> r_VEOG = %4.2f, r_HEOG = %4.2f, r_REOG = %4.2f, r_EKG = %4.2f\t'], ic.acc(ii, [1 2 4 3]));
	set(ff.tx, 'String', str);

	% % % the text for correlation 2
	str = [];
	for k = 1 : 9
		str = strcat(str, '  R', num2str(k, '%1d'), '=', num2str(abs(round(ic.lcc(ii, k)*1000)), '%03d'));
	end
	set(ff.tx2, 'String', str);
	clear str;

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
	pZ = plt_topoPlot2Data(ic.Z(:, ii), loc2d, fp.tpRangeXY, fp.tpTmpPtsNum);
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
	if ~isfield(ic, 'fM')
		ic.fM = ceil(ic.xfs(end)/10)*10;
	end
	set(ff.pwAX, 'XLim', [0 ic.fM], 'XTick', [2 4 10 30 50 ic.fM], 'XTickLabel', num2str([2 4 10 30 50 ic.fM]'));

	% % % trl plot
% 	set(ff.tr, 'XData', ic.trv(:, ii), 'YData', 1 : ic.nm);
% 	set(ff.trAX, 'XLim', [0 1], 'XTIck', [0 0.5 1], 'YLim', [0 ic.nm + 1]);
	drawnow;
end %end of function
