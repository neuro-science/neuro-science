function errCode = meg_raw_prepro_ica2clean(p, f, sti, rsp, is)
	% % % updated 26/09/2022 by wp for ica of clean data
	% % % The first datasets were for TP4
	%% 1. initialize
	% % % 1.1 prepare outputs	
	errCode = 0;

	% % % 1.2 check inputs	
	if nargin < 3
		errCode = sprintf( ...
			'Three parameters are needed: para, files and id');
		return;
	elseif ~exist(f.aus{is}, 'file')
		errCode = sprintf( ...
			'Input file %s was not found, exiting... \n', f.ein{is});
		return;
% 	elseif exist(f.skp{is}, 'file')
% 		errCode = sprintf( ...
% 			'The data %s is to skip, exiting... \n', f.idx{is});
% 		return;
	elseif exist(f.on{is}, 'file')
		errCode = sprintf( ...
			'The data %s is being processed, exiting... \n', f.idx{is});
		return;
	elseif exist(f.off{is}, 'file')
		errCode = sprintf( ...
			'The data %s was processed before, exiting... \n', f.idx{is});
		return;
	elseif ~exist(f.ica{is}, 'file')
		errCode = sprintf( ...
			'ICA file %s exist not, exiting... \n', f.ica{is});
		return;
	end
	
	%% 2. prepare para etc.
	% % % 2.1 place holder and welcome	
	save(f.on{is}, 'is');
	plt_myPrintLine([]);
	ttt = clock;
	fprintf('Process started on %s @%02d:%02d:%02.0f ...\n', ...
		f.idx{is}, ttt(4:6));

	% % % 2.2 get ica data	
	if strcmp(f.mic.oc{is}.fname, f.ica{is})
		v = load(f.ica{is});
		if strcmp(v.c.fname, f.mic.oc{is}.fname) && ...
				f.mic.oc{is}.n == v.c.nICs
			MM = v.c.A(:, ~f.mic.oc{is}.flag) * v.c.W(~f.mic.oc{is}.flag, :);
			clear v;
		else
			errCode = sprintf( ...
				'ICA data %s match not, exiting... \n', f.ica{is});
			return;
		end
	else
		errCode = sprintf( ...
			'ICA file %s match not, exiting... \n', f.ica{is});
		return;
	end
	
	% % % 2.3 preapare segmentation
	tmp = int32(p.pre.rawSampleRate * p.pro.dataInSecs(1)) : int32(p.pre.rawSampleRate * p.pro.dataInSecs(2));
	d0_id = tmp(1 : 3 : end)';
	npts = numel(d0_id);
% 	p.npts_b = length(b0_id);
	tmp = ismember(sti{is}(2, :), p.pro.onTrigger);
	t0 = int32(sti{is}(1, tmp));
	did = bsxfun(@plus, t0,	d0_id);
	clear d0_id t0 tmp;

	%% 3. do the data
	% % % 3.1 load data	
	clear meg;
	load(f.aus{is}, 'meg');
	ttt = clock;
	nchs = size(meg, 2);
	fprintf('Data loaded @%02d:%02d:%02.0f!\n', ttt(4:6));
	
	% % % 3.2 apply data	
% 	meg = (MM * v.meg')';
	meg = meg * MM';
	clear v;
	meg = reshape(meg(did, :), npts, [], nchs);
	
	% % % 3.3 save data	
	save(f.cln{is}, 'meg');
	clear meg;
	ttt = clock;
	fprintf('data saved @%02d:%02d:%02.0f!\n', ttt(4:6));
	plt_myPrintLine([]);
	save(f.off{is}, 'is');
	
end

