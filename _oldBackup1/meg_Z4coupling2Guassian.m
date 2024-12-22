function  out = meg_Z4coupling2Guassian (in, nTrls, nTapers, frqDim, trlDim, sFlag)

% % % updated 15/08/17 by wp: updated to new names
% % % updated 29/03/16 by wp: updated v by multiplying 2

% % % updated 25/11/14 by wp:
% % %		-	offer an option to set same channel coherence as 0
% % %		-	input check reformulated
% % % re-written 18/05/14 by wp
% % % for more flexible data

	%% 1. check headers
	sz = size(in);
	if nargin < 3 || isempty(nTapers)
		nTapers = 1;
	end
	
	if nargin < 4 || isempty(frqDim)
		frqDim = length(sz);
	end
	nFs = sz(frqDim);
	
	if nargin < 5 || isempty(trlDim)
		trlDim = length(sz);
	end
		
	if nargin < 6 || isempty(sFlag)
		sFlag = 0;
	end
		
	if sFlag
		if sz(1) == sz(2)
			in = reshape(in, [sz(1)*sz(2), prod(sz(3:end))]);
			in(1 : 1+sz(1) : end, :) = 0; % sent diagnal 0
			in = reshape(in, sz);
		else
			fprintf('First two dimensions mismatch, diagnal zero failed!\n');
		end
	end
	% % % coherence value > 1, could be cauculating inacurately 
	% % % because of limitation of single/double precisions.
	idx = find(in >= 1);
	if ~isempty(idx)
		r1 = 100 * length(idx) / numel(in);
		r2 = (max(in(idx)) - 1) * 100;
		if r2 < 1
			in(idx) = 1 - 1e-11;
			fprintf('!!! %5.1f%% of values are larger than 1, maximum %5.1f%%. !!!\n', r1, r2);
		else
			error('Incredibly coherence values, >1 by %6.1f%%.', r2);
		end
	end
	
	id0 = find(in < -1, 1);
	if ~isempty(id0)
		error('Incredibly coherence values, < -1 !');
	end
	
	%% 2. check whether tapers are not scalar
	if numel(nTapers) > 1
		nFts = length(nTapers);
		if nFts ~= nFs
			error('wrong size of frequencies!');
		else
			sz0 = ones(size(sz));
			sz0(frqDim) = sz(frqDim);
			sz1 = sz;
			sz1(frqDim) = 1;
			v = reshape(nTapers, sz0);
			vf = repmat(v, sz1);
			clear v sz0 sz1;
		end
	else
		vf = nTapers;
	end
	if numel(nTrls) > 1
		sz0 = ones(size(sz));
		sz0(trlDim) = sz(trlDim);
		sz1 = sz;
		sz1(trlDim) = 1;
		v = reshape(nTrls, sz0);
		vr = repmat(v, sz1);
		clear v sz0 sz1;
	else
		vr = nTrls;
	end
	v = 2 * vf .* vr; %%%add *2 on 29/03/2017
	
	%% 3. calculation
	out = 23 / 20 * (sqrt(log(1 - in .^ 2) .* (2 - v)) - 23 / 20);
	clear in v;
end

