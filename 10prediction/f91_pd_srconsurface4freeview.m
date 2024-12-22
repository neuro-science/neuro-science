function [dl, dr, ana] = f91_pd_srconsurface4freeview(data, ana, fname, para)
% % % updated 16/04/18 by wp: now we use oct3 with reduction (1154) nodes
% % % updated 20/12/17 by wp: now we use smaller number of grids
% % % written 15/08/17 by wp: generate activation file for freeview in pd project

	%% 0. decide where to put the data
	if nargin < 1
		fprintf('Input Error!\nWe need at least 1 inputs: SourceData!\n');
		return;
	end
		
	if ismac
		fpath = '/Users/wang/Documents/FreeSurfer/fsaverage5/';
		MNE_path = '/Applications/MNE_274_3420/share/matlab/';
	elseif isunix
		fpath = '/home/pwang/0FS/fsaverage5/';
		MNE_path = '/home/pwang/8TL/MNE274/share/matlab/';
	else
		error('We support only mac and linux');
	end
	
	%% 1. check inputs and prepare
	if nargin < 4 || isempty(para)
		para = 2;
	end
	if nargin < 3 || isempty(fname)
		fname = [fpath, 'surf/_tmp'];
	elseif isnumeric(fname)
		fname = [fpath, 'surf/_tmp', num2str(fname, '%02d')];
	else
		fname = [fpath, 'surf/_', fname];
	end
	lname = [fname, '_l.act'];
	rname = [fname, '_r.act'];
	
	%% 2. put data into source locations
	switch para
		case 3
			N1 = numel(ana(1).id3);
			N2 = numel(ana(2).id3);
			if max(abs(size(data) - [N1+N2 1])) > 0.1
				fprintf('Input Error!\nI support only 1154x1 voxels!\n');
				return;
			else
				ld0 = data(1 : N1);
				rd0 = data(N1 + 1 : N1 + N2);
			end
			dl = plt_srconsurface4freeview_NO_gaussian(lname, ld0, ana(1).id3, ana(1).vc, ana(1).ntri);	
			dr = plt_srconsurface4freeview_NO_gaussian(rname, rd0, ana(2).id3, ana(2).vc, ana(2).ntri);	
		case 2
			N1 = numel(ana(1).id2);
			N2 = numel(ana(2).id2);
			if max(abs(size(data) - [N1+N2 1])) > 0.1
				fprintf('Input Error!\nI support only 306x1 voxels!\n');
				return;
			else
				ld0 = data(1 : N1);
				rd0 = data(N1 + 1 : N1 + N2);
			end
			dl = plt_srconsurface4freeview_NO_gaussian(lname, ld0, ana(1).id2, ana(1).vc, ana(1).ntri);	
			dr = plt_srconsurface4freeview_NO_gaussian(rname, rd0, ana(2).id2, ana(2).vc, ana(2).ntri);	
		otherwise
			error('Only oct2 and oct3 are supported!\n');
		end

	
end
