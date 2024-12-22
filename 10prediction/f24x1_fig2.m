function f24x1_fig2

	%% 1. headers
	% % % 1.1 clean up the space
	clear all; close all; clc; tic; dbstop error;

	% % % 1.2 path settings
	if ismac
		fpath = '/Users/wang/Documents/6PD/';
		expath = '/Users/wang/ukebox/5Drafts/60PD_Prediction/3Source/';
	elseif isunix
		fpath = '~/6PD/';
		expath = '/mnt/homes/home020/pwang/4TP/5log/tmpFigs/';
	else
		fpath = 'D:\Data\6PD\';
		expath = 'C:\Users\peng_\ukebox\5Drafts\60PD_Prediction\3Source\';
	end
	
	% % % 1.3 load data
	load(fullfile(fpath, '6grp', 'on.mat'));
	load([fpath, '3sum/37neu/Plots/bar2scatter.mat'], 'bhv');	%[V/A; P/C; SUB, CR/RT]
	
	% % % 1.4 parameters of data
	sub_types = p.bhv.control_type;
	n1 = sum(sub_types == 1);
	n2 = sum(sub_types == 2);
	n = n1 + n2;
	
	% % % 1.5 parameters of plotting
	X0 = [3 6]; %X axis for two condition groups
	D0 = 1;	%the full deviation for conditions
	D1 = 0.1;	%the half deviation for cases
	markerSizes = [100, 100];	%Markersize
	markerColors = lines(n);
	markerStyles = {'o', '+', '^', '*', 'd', 's', 'p', 'h', 'x', '.'};
	markerThickness = 2;

	% % % 1.6 export to txt
% 	tmp = reshape(permute(bhv(:, :, :, 2), [3 1 2]), [], 4);
% 	fid = fopen([fpath, '3sum/37neu/Plots/rt.txt'], 'w');
% 	for k = 1 : p.nsb
% 		fprintf(fid, '%12.6f,', tmp(k, :));
% 		fprintf(fid, '\n');
% 	end
% 	fclose(fid);

	%% 2. plot 2A
	figure(101); 
	clf; hold on;
	set(gcf, 'color', [1 1 1]);

	% Plot each condition groups
	for i = 1:2	%P/C
		for j = 1:2 %V/A
			data = squeeze(bhv(j, i, :, 1));
			d1 = data(sub_types == 1);
			d2 = data(sub_types == 2);
			m = mean(data);
% 			m1 = mean(d1);
% 			m2 = mean(d2);
			se = std(data)/sqrt(n);

			% Plot males (cross) and females (circle)
			scatter(X0(i) + (j-1.5) * D0 + randn(1, n1) * D1, d1, ...
				markerSizes(1), markerColors(j, :), markerStyles{1}, 'LineWidth', markerThickness);
			scatter(X0(i) + (j-1.5) * D0 + randn(1, n2) * D1, d2, ...
				markerSizes(2), markerColors(j, :), markerStyles{2}, 'LineWidth', markerThickness);

% 			scatter(X0(i) + (j-1.5) * D0 + randn(1, n1) * D1, d1, ...
% 				markerSizes(1), markerColors(j, :), markerStyles{1}, 'filled');
% 			scatter(X0(i) + (j-1.5) * D0 + randn(1, n2) * D1, d2, ...
% 				markerSizes(2), markerColors(j, :), markerStyles{2}, 'filled');

			% Means with error bars (SEM)
			plot([X0(i) + (j-1.5) * D0 - D1, X0(i) + (j-1.5) * D0 + D1], [m, m], ...
				'-', 'color', markerColors(j, :), 'LineWidth', 3);
% 			plot([X0(i) + (j-1.5) * D0 - D1, X0(i) + (j-1.5) * D0 + D1], [m1, m1], ...
% 				'-.', 'color', markerColors(j, :), 'LineWidth', 1);
% 			plot([X0(i) + (j-1.5) * D0 - D1, X0(i) + (j-1.5) * D0 + D1], [m2, m2], ...
% 				'--', 'color', markerColors(j, :), 'LineWidth', 1);
			h = errorbar(X0(i) + (j-1.5) * D0, m, se, 'LineWidth', 1.5);
			set(h, 'color', markerColors(j, :));
		end
	end
	
	% Add significance markers manually
	text(4.5, 1.05, '*', 'HorizontalAlignment', 'center', 'fontsize', 30);
	plot([3, 6], [1.03, 1.03], 'k-', 'LineWidth', 1); % Line for significance

	% Aesthetics
	set(gca, 'fontsize', 36, ...
		'XTick', X0, 'XTickLabel', {'Prediction', 'Memory'}, ...
		'yLim', [0.55 1.1], 'YTick', [0.6 0.8 1], 'YTickLabel', [60 80 100]);
	ylabel('Corrected Accuracy (%)');
	export_fig(fullfile(expath, 'scatter_ca.eps'));

	
	%% 3. plot 2B
	figure(102); 
	clf; hold on;
	set(gcf, 'color', [1 1 1]);

	% Plot each condition groups
	for i = 1:2	%P/C
		for j = 1:2 %V/A
			data = squeeze(bhv(j, i, :, 2));
			d1 = data(sub_types == 1);
			d2 = data(sub_types == 2);
			m = mean(data);
% 			m1 = mean(d1);
% 			m2 = mean(d2);
			se = std(data)/sqrt(n);

			% Plot males (cross) and females (circle)
			scatter(X0(i) + (j-1.5) * D0 + randn(1, n1) * D1, d1, ...
				markerSizes(1), markerColors(j, :), markerStyles{1}, 'LineWidth', markerThickness);
			scatter(X0(i) + (j-1.5) * D0 + randn(1, n2) * D1, d2, ...
				markerSizes(2), markerColors(j, :), markerStyles{2}, 'LineWidth', markerThickness);

% 			scatter(X0(i) + (j-1.5) * D0 + randn(1, n1) * D1, d1, ...
% 				markerSizes(1), markerColors(j, :), markerStyles{1}, 'filled');
% 			scatter(X0(i) + (j-1.5) * D0 + randn(1, n2) * D1, d2, ...
% 				markerSizes(2), markerColors(j, :), markerStyles{2}, 'filled');

			% Means with error bars (SEM)
			plot([X0(i) + (j-1.5) * D0 - D1, X0(i) + (j-1.5) * D0 + D1], [m, m], ...
				'-', 'color', markerColors(j, :), 'LineWidth', 3);
% 			plot([X0(i) + (j-1.5) * D0 - D1, X0(i) + (j-1.5) * D0 + D1], [m1, m1], ...
% 				'-.', 'color', markerColors(j, :), 'LineWidth', 1);
% 			plot([X0(i) + (j-1.5) * D0 - D1, X0(i) + (j-1.5) * D0 + D1], [m2, m2], ...
% 				'--', 'color', markerColors(j, :), 'LineWidth', 1);
			h = errorbar(X0(i) + (j-1.5) * D0, m, se, 'LineWidth', 1.5);
			set(h, 'color', markerColors(j, :));
		end
	end
	

	% Add significance markers manually
	text(4.5, 2650, '*', 'HorizontalAlignment', 'center', 'fontsize', 30);
	plot([3, 6], [2575, 2575], 'k-', 'LineWidth', 1); % Line for significance

	text(3, 2450, '**', 'HorizontalAlignment', 'center', 'fontsize', 30);
	plot([2.5, 3.5], [2400, 2400], 'k-', 'LineWidth', 1); % Line for significance
	
	text(6, 2450, '***', 'HorizontalAlignment', 'center', 'fontsize', 30);
	plot([5.5, 6.5], [2400, 2400], 'k-', 'LineWidth', 1); % Line for significance

	% Aesthetics
	set(gca, 'fontsize', 36, ...
		'XTick', X0, 'XTickLabel', {'Prediction', 'Memory'}, ...
		'yLim', [0 2750], 'YTick', [0 1000 2000], 'YTickLabel', [0 1000 2000]);
	ylabel('Response Time (ms)');
	export_fig(fullfile(expath, 'scatter_rt.eps'));

	
	
	
	
	
	
	
	
	
end
