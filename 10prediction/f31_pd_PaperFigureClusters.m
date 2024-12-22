function f31_pd_PaperFigureClusters
	% % % all figures in the paper related to clusters

	%% 0. prepare environment
	% % % 0.1	clean up the space
	clear all; close all; clc; tic; dbstop error;
	
	% % % read para data and reset
	if ismac
		thePath = '/Users/wang/ukebox/5Drafts/60PD_Prediction/';
		fpath = '/Users/wang/Documents/6PD/';
	else
		thePath = '/Users/wang/ukebox/5Drafts/60PD_Prediction/';
		fpath = 'C:\Users\pwang\Documents\Data\';
	end
	load(fullfile(fpath, '6grp', 'on.mat'));
	
	
	%% 3. Figure 3 - evoked power 
% 	% % % 3.1 load data
% 	load(fullfile(thePath, '3Source', 'Fig3.mat'), 'f3');	
% 	
% 	% % % 3.2 spatial	
% 	% % % 3.2.1 export
% 	f91_pd_srconsurface4freeview(f3.spt1/1500, p.ana, 'ep_sp');
% 	% % % 3.2.2 run s08_freeview4PaperFigures
% 	% % % 3.2.3 merge
% 	f605_merge_plots('/Users/wang/Documents/6PD/4log/42figPaper/ccc/');
% 	
% 	% % % 3.3 tf - take care about the data aspectratio and axis size
% 	close all;
% 	cfg = [];
% 	cfg.fid = figure('position', [500 500 1000 1000]);
% 	set(gcf);
% 	cfg.szInterpolate = 10;
% 	cfg.cLim = 1000;
% 	cfg.aid = axes('fontsize', 36, 'parent', cfg.fid);
% 	cfg.xPos = [1 7 29];
% 	cfg.xAxisLabel = num2str([0 150 700]');
% 	cfg.th = 0;%.2;
% 	cfg = plt_tfPlot2d(f3.tft, cfg);
% 	set(cfg.aid, 'ytick', 3:4:19, 'yticklabel', [8, 16, 32, 64, 128], ...
% 		'FontName', 'helvetica', 'FontSize', 24, ...
% 		'position', [0.1, 0.1, 0.5, 0.3221]);
% 	xlabel('Time (ms)');
% 
% 	xt = linspace(0, 700, cfg.szInterpolate * (size(f3.tft, 2) - 1) + 1);
% 	yt = interpn(0:25:700, sum(f3.tft, 1), xt, 'pchip');
% 	at = axes(cfg.fid, 'position', [0.1, 0.4421, 0.5, 0.1]);
% 	h = plot(at, xt, yt, 'k', 'linewidth', 5);
% 	set(at, 'xlim', [0 700], 'xtick', [0 150 700], 'xticklabel', '', 'ylim', [0 15000], 'ytick', [0 15000], 'FontName', 'helvetica', 'FontSize', 24);
% 	xlabel('');
% 
% 			
% 	xf = linspace(2.5, 7, (size(f3.tft, 1) - 1) * 10 + 1);
% 	yf = interpn(linspace(2.5, 7, size(f3.tft, 1)), sum(f3.tft, 2), xf, 'pchip');
% 	af = axes(cfg.fid, 'position', [0.62, 0.1, 0.1, 0.3221]);
% 	h = plot(af, yf, xf, 'k', 'linewidth', 5);
% 	set(af, 'xlim', [0 15000], 'ytick', 3:7, 'yticklabel', '', 'ylim', [2.5 7], 'xtick', [0 15000], 'FontName', 'helvetica', 'FontSize', 24);
% 			
% 	export_fig(fullfile(thePath,  '4Inter', 'Figure3', 'f3_ep_tf.eps'));
% 
% 	% % % 3.4 color bar
% 	clf;
% 	y = plt_scale2color(linspace(-1, 1, 101), 0, 1);
% 	y = permute(repmat(y, [1 1 20]), [1 3 2]);
% 	image(y)
% 	axis equal;
% 	export_fig(fullfile(thePath,  '4Inter', 'Figure3', 'f3_bar.eps'));
% 	image(y(1:51, :, :));
% 	axis equal;
% 	export_fig(fullfile(thePath,  '4Inter', 'Figure3', 'f3_bar-.eps'));
% 	image(y(51:101, :, :));
% 	axis equal;
% 	export_fig(fullfile(thePath,  '4Inter', 'Figure3', 'f3_bar+.eps'));


	%% 4. Figure 4 - power cluster for anova
% 	% % % 4.1 load data
% 	load(fullfile(thePath, '3Source', 'Fig4.mat'), 'f4');	
% 	
% 	% % % 4.2 spatial	
% 	% % % 4.2.1 export
% 	for ic = 1 : 4
% 		f91_pd_srconsurface4freeview(f4.spf(:, ic)/f4.M1(ic), p.ana, [f4.tag{ic}, '_sp']);
% 	end
% 	% % % 4.2.2 run s08_freeview4PaperFigures
% 	% % % 4.2.3 merge
% 	f605_merge_plots('/Users/wang/Documents/6PD/4log/42figPaper/ccc/');
% 	
% 	close all;
% 	N = 10;
% 	Dx = [0 8000];
% 	Dy = [0 6000];
% 	for ic = 1 : 4
% 		cfg = [];
% 		cfg.fid = figure('position', [500 500 1000 1000]);
% 		cfg.szInterpolate = N;
% 		cfg.cLim = f4.M2(ic);
% 		cfg.aid = axes('fontsize', 36, 'parent', cfg.fid);
% 		cfg.xPos = [1 6];
% 		cfg.xAxisLabel = num2str([100 600]');
% 		cfg = plt_tfPlot2d(f4.tff(:, :, ic), cfg);
% 		set(cfg.aid, 'DataAspectRatio', [1 2 1], 'ytick', 3:4:19, 'yticklabel', [8 16 32 64 128], ...
% 			'FontName', 'helvetica', 'FontSize', 24, 'position', [0.1, 0.1, 0.363, 0.6442]);
% 		xlabel('Time (ms)');
% 		
% 		xt = linspace(p.pro.ie1.xts1(1), p.pro.ie1.xts1(end), (size(f4.tff(:, :, ic), 2) - 1) * N + 1);
% 		yt = interpn(p.pro.ie1.xts1, sum(f4.tff(:, :, ic), 1), xt, 'pchip');
% 		at = axes(cfg.fid, 'position', [0.1, 0.7642, 0.363, 0.1]);
% 		h = plot(at, xt, yt, 'k', 'linewidth', 5);
% 		set(at, 'xlim', [100 600], 'xtick', [100 600], 'xticklabel', '', 'ylim', Dx, 'ytick', Dx, 'FontName', 'helvetica', 'FontSize', 24);
% 		ylabel('\itF\rm-value');
% 		
% 		xf = linspace(2.5, 7, (size(f4.tff(:, :, ic), 1) - 1) * N + 1);
% 		yf = interpn(linspace(2.5, 7, 19), sum(f4.tff(:, :, ic), 2), xf, 'pchip');
% 		af = axes(cfg.fid, 'position', [0.483, 0.1, 0.1, 0.6442]);
% 		h = plot(af, yf, xf, 'k', 'linewidth', 5);
% 		set(gca, 'xlim', Dy, 'ytick', 3:7, 'yticklabel', '', 'ylim', [2.5 7], 'xtick', Dy, 'FontName', 'helvetica', 'FontSize', 24);
% 		xlabel('\itF\rm-value');
% 		export_fig(fullfile(thePath,  '4Inter', 'Figure4', ['f4', f4.tag{ic}, '_tf.eps']));
% 	end
	
	
	%% 6. Figure 5 - imc cluster
% 	% % % 4.1 load data
% 	load(fullfile(thePath, '3Source', 'Fig5.mat'), 'f5');	
% 	
% 	% % % 4.2 spatial	
% 	% % % 4.2.1 export
% 	for ic = 1 : 6
% 		f91_pd_srconsurface4freeview(f5.sp(:, ic)/f5.M1(ic), p.ana, [f5.tag{ic}, '_sp']);
% 	end
% 	% % % 4.2.2 run s08_freeview4PaperFigures
% 	% % % 4.2.3 merge
% % 	f605_merge_plots('/Users/wang/Documents/6PD/4log/42figPaper/ccc/');
% 	
% 	close all;
% 	N = 10;
% 	for ic = 1 : 6
% 		cfg = [];
% 		cfg.fid = figure('position', [500 500 1000 1000]);
% 		cfg.szInterpolate = N;
% 		cfg.cLim = f5.M2(ic);
% 		cfg.aid = axes('fontsize', 36, 'parent', cfg.fid);
% 		cfg.xPos = [1 6];
% 		cfg.xAxisLabel = num2str([100 600]');
% 		cfg = plt_tfPlot2d(f5.tf(:, :, ic), cfg);
% 		set(cfg.aid, 'DataAspectRatio', [1 2 1], 'ytick', 3:4:19, 'yticklabel', [8 16 32 64 128], ...
% 			'FontName', 'helvetica', 'FontSize', 24, 'position', [0.1, 0.1, 0.363, 0.6442]);
% 		xlabel('Time (ms)');
% 		
% 		xt = linspace(p.pro.ie1.xts1(1), p.pro.ie1.xts1(end), (size(f5.tf(:, :, ic), 2) - 1) * N + 1);
% 		yt = interpn(p.pro.ie1.xts1, sum(f5.tf(:, :, ic), 1), xt, 'pchip');
% 		at = axes(cfg.fid, 'position', [0.1, 0.7642, 0.363, 0.1]);
% 		h = plot(at, xt, yt, 'k', 'linewidth', 5);
% 		set(at, 'xlim', [100 600], 'xtick', [100 600], 'xticklabel', '', 'ylim', f5.dx(ic, :), 'ytick', f5.dx(ic, :), 'FontName', 'helvetica', 'FontSize', 24);
% 		ylabel('\itF\rm-value');
% 		
% 		xf = linspace(2.5, 7, (size(f5.tf(:, :, ic), 1) - 1) * N + 1);
% 		yf = interpn(linspace(2.5, 7, 19), sum(f5.tf(:, :, ic), 2), xf, 'pchip');
% 		af = axes(cfg.fid, 'position', [0.483, 0.1, 0.1, 0.6442]);
% 		h = plot(af, yf, xf, 'k', 'linewidth', 5);
% 		set(gca, 'xlim', f5.dy(ic, :), 'ytick', 3:7, 'yticklabel', '', 'ylim', [2.5 7], 'xtick', f5.dy(ic, :), 'FontName', 'helvetica', 'FontSize', 24);
% 		xlabel('\itF\rm-value');
% 		export_fig(fullfile(thePath,  '4Inter', 'Figure6', ['f5', f5.tag{ic}, '_tf.eps']));
% 	end


	
	%% 9. Figure 9 - history and future
% 	% % % 9.1 prepare data and paras	
% 	load(fullfile(thePath, '3Source', 'Fig7.mat'), 'f7');	
% 	
% 	% % % 9.2 spatial	
% 	% % % 9.2.1 export
% 	for ic = 1 : 2
% 		f91_pd_srconsurface4freeview(f7.sp(:, ic)/f7.M1, p.ana, [f7.str{ic}, '_sp']);
% 	end
% 	% % % 9.2.2 run s08_freeview4PaperFigures
% 	% % % 9.2.3 merge
% % 	f605_merge_plots('/Users/wang/Documents/6PD/4log/45figPaper/ccc/');
% 	
% 	% % % 9.3 tf
% 	close all;
% 	N = 10;
% 	for ic = 1 : 2
% 		cfg = [];
% 		cfg.fid = figure('position', [500 500 1000 1000]);
% 		cfg.szInterpolate = N;
% 		cfg.cLim = f7.M2;
% 		cfg.aid = axes('fontsize', 36, 'parent', cfg.fid);
% 		cfg.xPos = [1 6];
% 		cfg.xAxisLabel = num2str([100 600]');
% 		cfg = plt_tfPlot2d(f7.tf(:, :, ic), cfg);
% 		set(cfg.aid, 'DataAspectRatio', [1 2 1], 'ytick', 3:4:19, 'yticklabel', [8 16 32 64 128], ...
% 			'FontName', 'helvetica', 'FontSize', 24, 'position', [0.1, 0.1, 0.363, 0.6442]);
% 		xlabel('Time (ms)');
% 		
% 		xt = linspace(p.pro.ie1.xts1(1), p.pro.ie1.xts1(end), (size(f7.tf(:, :, ic), 2) - 1) * N + 1);
% 		yt = interpn(p.pro.ie1.xts1, sum(f7.tf(:, :, ic), 1), xt, 'pchip');
% 		at = axes(cfg.fid, 'position', [0.1, 0.7842, 0.363, 0.1]);
% 		h = plot(at, xt, yt, 'k', 'linewidth', 5);
% 		set(at, 'xlim', [100 600], 'xtick', [100 600], 'xticklabel', '', 'ylim', f7.dx, 'ytick', f7.dx, 'FontName', 'helvetica', 'FontSize', 24);
% 		ylabel('\itF\rm-value');
% 		
% 		xf = linspace(2.5, 7, (size(f7.tf(:, :, ic), 1) - 1) * N + 1);
% 		yf = interpn(linspace(2.5, 7, 19), sum(f7.tf(:, :, ic), 2), xf, 'pchip');
% 		af = axes(cfg.fid, 'position', [0.503, 0.1, 0.1, 0.6442]);
% 		h = plot(af, yf, xf, 'k', 'linewidth', 5);
% 		set(gca, 'xlim', f7.dy, 'ytick', 3:7, 'yticklabel', '', 'ylim', [2.5 7], 'xtick', f7.dy, 'FontName', 'helvetica', 'FontSize', 24);
% 		xlabel('\itF\rm-value');
% 		export_fig(fullfile(thePath,  '4Inter', 'Figure9', ['f7_', f7.str{ic}, '_tf.eps']));
% 	end
% 	
% 	% % % 9.4 modulation
% 	LW = 5;
% 	for ic = 1 : 2
% 		close all;
% 		h1 = bar(f7.m(:, :, ic), 'LineWidth', LW);
% 		set(h1(1), 'facecolor', [1 1 1]*0.2);
% 		set(h1(2), 'facecolor', [1 1 1]*0.8);
% 		hold on;
% 		h2 = errorbar(f7.m(:, :, ic), f7.e(:, :, ic), 'LineStyle', 'none', 'color', [0 0 0], 'linewidth', LW);
% 		drawnow;
% 		set(gca, 'fontsize', 36, 'box', 'off', 'LineWidth', LW, 'ytick', [-0.02 0 0.02 0.04], 'yticklabel', {'-2', '0', '2', '4'}, 'ylim', [-0.02 0.04], 'xticklabel', {'Sampling Period', 'Predicting Period'});
% 		fname = fullfile(thePath,  '4Inter', 'Figure7', ['f7_bar_', f7.str{ic}, '.eps']);
% 		export_fig(fname);	
% 	end	

	%% 11. Figure S3 - response power
% 	% % % 11.1 prepare data and paras	
% 	load(fullfile(thePath, '3Source', 'FigS3.mat'), 'f11');	
% 	
% 	% % % 11.2 spatial	
% 	% % % 11.2.1 export
% 	for ic = 1 : 3
% 		f91_pd_srconsurface4freeview(f11.sp(:, ic)/f11.M1, p.ana, ['rsp_ip', num2str(ic), '_sp']);
% 	end
% 	% % % 11.2.2 run s08_freeview4PaperFigures
% 	% % % 11.2.3 merge
% % 	f605_merge_plots('/Users/wang/Documents/6PD/4log/45figPaper/ccc/');
% 	
% 	% % % 11.3 tf
% 	close all;
% 	N = 10;
% 	for ic = 1 : 3
% 		cfg = [];
% 		cfg.fid = figure('position', [500 500 1000 1000]);
% 		cfg.szInterpolate = N;
% 		cfg.cLim = f11.M2;
% 		cfg.aid = axes('fontsize', 36, 'parent', cfg.fid);
% 		cfg.xPos = [1 6 13];
% 		cfg.xAxisLabel = num2str([100 600 1300]');
% 		cfg = plt_tfPlot2d(f11.tf(:, :, ic), cfg);
% 		set(cfg.aid, 'DataAspectRatio', [1 1 1], 'ytick', 3:4:19, 'yticklabel', [8 16 32 64 128], ...
% 			'FontName', 'helvetica', 'FontSize', 24, 'position', [0.1, 0.1, 0.4307, 0.6442]);
% 		xlabel('Time (ms)');
% 		
% 		xt = linspace(p.pro.ie1.xts2(1), p.pro.ie1.xts2(end), (size(f11.tf(:, :, ic), 2) - 1) * N + 1);
% 		yt = interpn(p.pro.ie1.xts2, sum(f11.tf(:, :, ic), 1), xt, 'pchip');
% 		at = axes(cfg.fid, 'position', [0.1, 0.7842, 0.4307, 0.1]);
% 		h = plot(at, xt, yt, 'k', 'linewidth', 5);
% 		set(at, 'xlim', [100 1300], 'xtick', [100 600 1300], 'xticklabel', '', 'ylim', f11.dx(ic, :), 'ytick', f11.dx(ic, :), 'FontName', 'helvetica', 'FontSize', 24);
% 		ylabel('\itt\rm-value');
% 		
% 		xf = linspace(2.5, 7, (size(f11.tf(:, :, ic), 1) - 1) * N + 1);
% 		yf = interpn(linspace(2.5, 7, 19), sum(f11.tf(:, :, ic), 2), xf, 'pchip');
% 		af = axes(cfg.fid, 'position', [0.5707, 0.1, 0.1, 0.6442]);
% 		h = plot(af, yf, xf, 'k', 'linewidth', 5);
% 		set(gca, 'xlim', f11.dy(ic, :), 'ytick', 3:7, 'yticklabel', '', 'ylim', [2.5 7], 'xtick', f11.dy(ic, :), 'FontName', 'helvetica', 'FontSize', 24);
% 		xlabel('\itt\rm-value');
% 		export_fig(fullfile(thePath,  '4Inter', 'Figure11', ['f11_', num2str(ic), '_tf.eps']));
% 	end
	
	%% 12. Figure 12 - response imc
% 	% % % 12.1 prepare data and paras	
% 	load(fullfile(thePath, '3Source', 'FigS4.mat'), 'f12');	
% 	
% 	% % % 12.2 spatial	
% 	% % % 12.2.1 export
% 	for ic = 1 : 2
% 		f91_pd_srconsurface4freeview(f12.sp(:, ic)/f12.M1(ic), p.ana, ['rsp_ic', num2str(ic), '_sp']);
% 	end
% 	% % % 12.2.2 run s08_freeview4PaperFigures
% 	% % % 12.2.3 merge
% 	f605_merge_plots('/Users/wang/Documents/6PD/4log/45figPaper/ccc/');
% 	
% 	% % % 12.3 tf
% 	close all;
% 	N = 10;
% 	for ic = 1 : 2
% 		cfg = [];
% 		cfg.fid = figure('position', [500 500 1000 1000]);
% 		cfg.szInterpolate = N;
% 		cfg.cLim = f12.M2(ic);
% 		cfg.aid = axes('fontsize', 36, 'parent', cfg.fid);
% 		cfg.xPos = [1 6 13];
% 		cfg.xAxisLabel = num2str([100 600 1300]');
% 		cfg = plt_tfPlot2d(f12.tf(:, :, ic), cfg);
% 		set(cfg.aid, 'DataAspectRatio', [1 1 1], 'ytick', 3:4:19, 'yticklabel', [8 16 32 64 128], ...
% 			'FontName', 'helvetica', 'FontSize', 24, 'position', [0.1, 0.1, 0.4307, 0.6442]);
% 		xlabel('Time (ms)');
% 		
% 		xt = linspace(p.pro.ie1.xts2(1), p.pro.ie1.xts2(end), (size(f12.tf(:, :, ic), 2) - 1) * N + 1);
% 		yt = interpn(p.pro.ie1.xts2, sum(f12.tf(:, :, ic), 1), xt, 'pchip');
% 		at = axes(cfg.fid, 'position', [0.1, 0.7842, 0.4307, 0.1]);
% 		h = plot(at, xt, yt, 'k', 'linewidth', 5);
% 		set(at, 'xlim', [100 1300], 'xtick', [100 600 1300], 'xticklabel', '', 'ylim', f12.dx(ic, :), 'ytick', f12.dx(ic, :), 'FontName', 'helvetica', 'FontSize', 24);
% 		ylabel('\itt\rm-value');
% 		
% 		xf = linspace(2.5, 7, (size(f12.tf(:, :, ic), 1) - 1) * N + 1);
% 		yf = interpn(linspace(2.5, 7, 19), sum(f12.tf(:, :, ic), 2), xf, 'pchip');
% 		af = axes(cfg.fid, 'position', [0.5707, 0.1, 0.1, 0.6442]);
% 		h = plot(af, yf, xf, 'k', 'linewidth', 5);
% 		set(gca, 'xlim', f12.dy(ic, :), 'ytick', 3:7, 'yticklabel', '', 'ylim', [2.5 7], 'xtick', f12.dy(ic, :), 'FontName', 'helvetica', 'FontSize', 24);
% 		xlabel('\itt\rm-value');
% 		export_fig(fullfile(thePath,  '4Inter', 'Figure12', ['f12_', num2str(ic), '_tf.eps']));
% 	end

