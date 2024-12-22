% % % This function is to remove 50Hz line artifact based on chronux functions
% % % written 27/09/2016 by wp for more general purpose
% % % Inputs:
% % %		data (points, chans) - the input noisy data
% % %		segLength - length for each segment
% % %		bad_idx - the indices for bad data (e.g. detected by visual inspection etc.)
function [data, out] = cmp58_rawmeg_rm_50Hz(data, segLength, Fs, f0, bad_idx, fPass, tapers, pp)

	%% 1.initialization & parameters
	tic;
	if nargin < 8 || isempty(pp)
		pp = 0.01;
	end
	if nargin < 7 || isempty(tapers)
		p.tapers = [3 5];
	else
		p.tapers = tapers;
	end
	if nargin < 6 || isempty(fPass)
		p.fpass = [3 5];
	else
		p.fpass = fPass;
	end
	if nargin < 5
		bad_idx = [];
	end
	if nargin < 4 || isempty(f0)
		f0 = 50;
	end
	if nargin < 3 || isempty(Fs)
		p.Fs = 1200;
	else
		p.Fs = Fs;
	end
	if nargin < 2 || isempty(segLength)
		segLength = p.Fs * 10;
	end
	if nargin < 1 || isempty(data)
		error('Input data cannot be empty!');
	end
	
	%% 2. prepare data
	[npts, nchs] = size(data); %[points, channels]
	idx = true(npts, 1);
	idx(bad_idx) = false;
	N = npts - length(bad_idx); %good points
	M = ceil(N / segLength);
	
	% % % save paras
	out = p;
	out.plt = 'n';
	out.pp = pp;
	out.f0 = f0;
	out.segLength = segLength;
	out.bad_idx = bad_idx;
	
	%% 3. work on data in loops (save memory)
	parfor ic = 1 : nchs
		% % % take data from one channel and fold into several segments
		d0 = data(idx, ic);
		d0(N + 1 : M * segLength) = 0;
		d0 = reshape(d0, [segLength, M]);
		% % % work on it
		d = rmlinesc(d0, p, pp, 'n', f0);
		% % % fill the real data
		d = d(:);
		d(N + 1 : M * segLength) = [];
		data(idx, ic) = d;
	end

	clear d d0;
	toc;
end