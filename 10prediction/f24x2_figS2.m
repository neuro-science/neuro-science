function f24x2_figS2

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
	% % %1-NV;2-NA;3-PV;4-PA;
	% % %1~4 sampling 5~8 predicting;

	
	% % % 1.5 parameters of data
	sub_types = p.bhv.control_type;
	n1 = sum(sub_types == 1);
	n2 = sum(sub_types == 2);
	n = n1 + n2;
	
	% % % 1.4 parameters of plotting
	X0 = [1, 3]; %X axis for two condition groups
	D1 = 0.3;	%the half deviation for cases
	markerSizes = [200, 200];	%Markersize
	markerColors = [143 195 31; 248 182 45]/255;
	markerStyles = {'o', '+', '^', '*', 'd', 's', 'p', 'h', 'x', '.'};
	markerThickness = 2;

	%% 2. plot S2A (power)
	figure(101); 
	clf; hold on;
	set(gcf, 'color', [1 1 1]);
	
	load([fpath, '3sum/37neu/data_pwr.mat'], 'c');	%[V/A; P/C; SUB, CR/RT]
	for i = 1:2
		% Plot each condition groups
		data = permute(mean(reshape(nanmean(c{i*2}.data(:, (1:4) + (i-1) * 4, :), 1), 2, 2, p.nsb), 1), [3 2 1]);
		d = data(:, 1) - data(:, 2);
		for j = 1:2
			d1 = d(sub_types == j);

			scatter(X0(i) + randn(1, numel(d1)) * D1, d1, ...
				markerSizes(i), markerColors(i, :), markerStyles{j}, 'LineWidth', markerThickness);

		end
	end
	plot([X0(1) - 1, X0(end) + 1], [0 0], 'k-');

	% Aesthetics
	set(gca, 'fontsize', 20, 'XTick', X0, 'XTickLabel', {'Sampling Segments', 'Prediction Segments'}, ...
		'yLim', [-5000 15000], 'YTick', [-5000 0 5000 10000 15000], 'YTickLabel', [-5000 0 5000 10000 15000]);
	ylabel('Arbitray unit');
	xlabel('Beta Power (Memory - Prediction)');
	export_fig(fullfile(expath, 'pwr_sub.eps'));
	clear c data d d1;
	
	%% 3. plot S2B (coupling theta)
	load([fpath, '3sum/37neu/data_imc.mat'], 'c');	%[V/A; P/C; SUB, CR/RT]

	figure(102); 
	clf; hold on;
	set(gcf, 'color', [1 1 1]);
	
	for i = 1:2
		% Plot each condition groups
		data = permute(nanmean(nanmean(c{i}.data(c{i}.tpe, :, :, :), 1), 3), [2 4 1 3]);
		d = data(:, 2) - data(:, 1);
		for j = 1:2
			d1 = d(sub_types == j);

			scatter(X0(i) + randn(1, numel(d1)) * D1, d1, ...
				markerSizes(i), markerColors(i, :), markerStyles{j}, 'LineWidth', markerThickness);

		end
	end
	plot([X0(1) - 1, X0(end) + 1], [0 0], 'k-');

	% Aesthetics
	set(gca, 'fontsize', 20, 'XTick', X0, 'XTickLabel', {'Sampling Segments', 'Prediction Segments'}, ...
		'yLim', [-0.005 0.04], 'YTick', [0 0.02 0.04], 'YTickLabel', [0 0.02 0.04]);
	ylabel('Arbitray unit');
	xlabel('Theta Coupling (Prediction - Memory)')
	export_fig(fullfile(expath, 'theta_sub.eps'));
	clear c data d d1;

	
	%% 4. plot S2C (coupling gamma)
	load([fpath, '3sum/37neu/data_imc.mat'], 'c');	%[V/A; P/C; SUB, CR/RT]

	figure(103); 
	clf; hold on;
	set(gcf, 'color', [1 1 1]);
	
	for i = 1:2
		% Plot each condition groups
		data = permute(nanmean(nanmean(c{i}.data(c{i}.tne, :, :, :), 1), 3), [2 4 1 3]);
		d = data(:, 1) - data(:, 2);
		for j = 1:2
			d1 = d(sub_types == j);

			scatter(X0(i) + randn(1, numel(d1)) * D1, d1, ...
				markerSizes(i), markerColors(i, :), markerStyles{j}, 'LineWidth', markerThickness);

		end
	end
	plot([X0(1) - 1, X0(end) + 1], [0 0], 'k-');

	% Aesthetics
	set(gca, 'fontsize', 20, 'XTick', X0, 'XTickLabel', {'Sampling Segments', 'Prediction Segments'}, ...
		'yLim', [-0.005 0.025], 'YTick', [0 0.01 0.02], 'YTickLabel', [0 0.01 0.02]);
	ylabel('Arbitray unit');
	xlabel('Gamma Coupling (Memory - Prediction)');
	export_fig(fullfile(expath, 'gamma_sub.eps'));
	clear c data d d1;

end
