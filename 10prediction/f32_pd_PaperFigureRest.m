function f32_pd_PaperFigureRest
	% % % all figures in the paper and related data

	%% 0. prepare environment
	% % % 0.1	clean up the space
	clear all; close all; clc; tic; dbstop error;
	
	% % % read para data and reset
	if ismac
		thePath = '/Users/wang/ukebox/5Drafts/60PD_Prediction/';
		fpath = '/Users/wang/Documents/6PD/';
	else
		thePath = 'C:\Users\pwang\ukebox\5Drafts\60PD_Prediction/';
		fpath = 'C:\Users\pwang\Documents\Data\';
	end
	load(fullfile(fpath, '6grp', 'on.mat'));
	
	%% 7. Figure 7 - phase coupling
% 	% % % 7.1 prepare data and paras	
% 	load(fullfile(thePath, '3Source', 'Fig6.mat'), 'f6');	
% 	N = 8;
% 	x = linspace(-1, 1, N+1);
% 	COLOR = [229 0 17; 0 153 67]/255;
% 	WIDTH = [2, 2];
% 	MARKER = {'+', 'x'};
% 	ALPHA = [0.1 0.1];
% 	
% 	% % % 7.2 cr/rt to phase
% 	close all;
% 	im = 1;
% 	D = [1e-3 1];
% 	for ib = 1 : 2
% 		for k1 = 1 : 2
% 			figure('PaperOrientation', 'landscape');
% 			clf;
% 			for k2 = 1 : 2
% 				h = shadedErrorBar(x, f6.m([1:N, 1], k2, k1, ib, im), f6.e([1:N, 1], k2, k1, ib, im));
% 				set(h.mainLine, 'Color', COLOR(k2, :), 'LineWidth', WIDTH(k2), 'Marker', MARKER{k2}, 'MarkerSize', 10);
% 				set(h.patch, 'FaceColor', COLOR(k2, :), 'FaceAlpha', ALPHA(k2));
% 			end
% 			set(gca, 'xlim', [-1 1], 'ylim', [-1 1]*D(ib), 'ytick', [-1 0 1]*D(ib), 'xtick', [-1 -0.5 0 0.5 1], 'xticklabel', {'-\pi' '-\pi/2' '0' '\pi/2' '\pi'}, 'fontsize', 60);
% 			fname = fullfile(thePath,  '4Inter', 'Figure6', ['f6_ic_ph_', num2str(ib), num2str(k1),'.pdf']);
% 			print(gcf, '-dpdf', '-fillpage', fname); 
% 		end
% 	end
% 	
% 	% % % 7.3 cr/rt to distance
% 	im = 2;
% 	D = [1e-3 1];
% 	close all;
% 	for ib = 1 : 2
% 		for k1 = 1 : 2
% 			figure;
% 			hold on;
% 			for k2 = 1 : 2
% 				px = polyfit((1:N)', f6.m(1:N, k2, k1, ib, im), 1);
% 				py = polyval(px, (1:N)');
% 				h = errorbar(f6.m(1:N, k2, k1, ib, im), f6.e(1:N, k2, k1, ib, im), '.', 'markersize', 60, 'color', COLOR(k2, :));
% 				plot((1:N)', py, '--', 'color', COLOR(k2, :), 'LineWidth', 2);
% 				[x1, x3(k2, k1, ib)] = corr((1:N)', f6.m(1:N, k2, k1, ib, im), 'Type','Spearman');
% 				str = sprintf('r = %4.2f', x1);
% 				text(4, max(f6.m(1:N, k2, k1, ib, im))*k2*2 - 3 * max(f6.m(1:N, k2, k1, ib, im)), str, 'fontsize', 24, 'color', COLOR(k2, :));
% 			end
% 			set(gca, 'xlim', [0 N+1], 'ylim', [-1 1]*D(ib), 'ytick', [-1 0 1]*D(ib), 'xtick', [1 8], 'xticklabel', {'min', 'max'}, 'fontsize', 60);
% 			fname = fullfile(thePath,  '4Inter', 'Figure6', ['f6_ic_pl_', num2str(ib), num2str(k1),'.eps']);
% 			export_fig(fname);
% 		end
% 	end
% 	
% 	% % % 7.4 significance with voxels
% 	for im = 1 : 2
% 		for it = 1 : 2
% 			figure;
% 			hold on;
% 			for ic = 1 : 2
% 				theM = 4;
% 				h = plot(f6.sp(:, ic, it, im, 3), f6.sp(:, ic, it, im, theM), '.', 'markersize', 20, 'color', COLOR(ic, :));
% 				px = polyfit(f6.sp(:, ic, it, im, 3), f6.sp(:, ic, it, im, theM), 1);
% 				py = polyval(px, 1 : 4000);
% 				plot((1 : 4000)', py, '--', 'color', COLOR(ic, :), 'LineWidth', 2);
% 				[x1, x2] = corrcoef(f6.sp(:, ic, it, im, 3), f6.sp(:, ic, it, im, theM));
% 				x3(ic, it, im) = x2(1, 2);
% 				if x3(ic, it, im) < 0.001
% 					text(2000, ic*200-300, ['{\itr} = ', num2str(x1(1, 2), '%4.2f'), '***'], 'fontsize', 60, 'color', COLOR(ic, :));
% 				elseif x3(ic, it, im) < 0.01
% 					text(2000, ic*200-300, ['{\itr} = ', num2str(x1(1, 2), '%4.2f'), '**'], 'fontsize', 60, 'color', COLOR(ic, :));
% 				elseif x3(ic, it, im) < 0.05
% 					text(2000, ic*200-300, ['{\itr} = ', num2str(x1(1, 2), '%4.2f'), '*'], 'fontsize', 60, 'color', COLOR(ic, :));
% 				else
% 					text(2000, ic*200-300, ['{\itr} = ', num2str(x1(1, 2), '%4.2f')], 'fontsize', 60, 'color', COLOR(ic, :));
% 				end
% 			end
% 			set(gca, 'xlim', [0 3500], 'ylim', [-1 1]*120, 'ytick', [-1 0 1]*100, 'xtick', [0 3000], 'fontsize', 60);
% 			fname = fullfile(thePath,  '4Inter', 'Figure6', ['f6_ic_pl_corr', num2str(im), num2str(it),'.eps']);
% 			export_fig(fname);
% 		end
% 	end		
	

	%% 10. Figure 8 - pac % later version figure 7 C, F
% 	load(fullfile(thePath, '3Source', 'Fig8.mat'), 'f8');	
% 	LW = 5; Y = [0.2 0.1];
% 	for ic = 1 : 2
% 		figure(ic + 100);clf;
% 		h1 = bar(f8.m1(:, :, ic), 'LineStyle', 'none');
% 		set(h1(1), 'facecolor', [1 1 1]*0.2);
% 		set(h1(2), 'facecolor', [1 1 1]*0.8);
% 		hold on;
% 		h2 = errorbar(f8.m1(:, :, ic), f8.e1(:, :, ic), 'LineStyle', 'none', 'color', [0 0 0], 'linewidth', LW);
% 		set(h2(1), 'color', [1 1 1]*0.2);
% 		set(h2(2), 'color', [1 1 1]*0.8);
% 		drawnow;
% 		set(gca, 'fontsize', 60, 'box', 'off', 'ytick', [0 Y(ic)], 'yticklabel', {'0', num2str(100*Y(ic))}, 'ylim', [0 Y(ic)], 'xticklabel', {'Sampling Period', 'Predicting Period'});
% 		export_fig(['/Users/wang/ukebox/5Drafts/60PD_Prediction/4Inter/Figure8/pac_pr_nx_', num2str(ic), '.eps']);
% 	end
		
end

