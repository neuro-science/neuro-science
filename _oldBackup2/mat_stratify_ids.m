function out = mat_stratify_ids (in, nid, NNN)
	% % % updated by wp @22/03/2018, stratify ids for N cons
	if nargin < 3
		NNN = 999;
	end
	if iscell(in)
		% % % get info
		N = numel(in);
		n = zeros(N, 1);
		for k = 1 : N
			in{k} = in{k}(:);
			n(k) = length(in{k});
		end
	elseif isvector(in)
		ns = unique(in);
		ns(ns > NNN) = [];	%remove extremely large ins
		N = numel(ns);
		if nargin > 1 && ~isempty(nid)
			in(nid) = N + NNN;
		end
		in2 = cell(N, 1);
		for k = 1 : N
			in2{k} = find(in == ns(k));
			n(k) = length(in2{k});
		end
		in = in2;
	else
		error('unsupported input, we support only cell and vector!');
	end
	
	
	% % % do change
	nn = min(n); %get the smallest number
	out = zeros(nn, N);
	for k = 1 : N
		out(:, k) = chooseHere(in{k}, nn);
	end
end

% function [out1, out2] = mat_stratify_ids (in1, in2)
% 	% % % first written by wp @18/12/2017, stratify ids for two cons
% 	in1 = in1(:);
% 	in2 = in2(:);
% 	n1 = length(in1);
% 	n2 = length(in2);
% 	if n1 <= n2
% 		out1 = in1;
% 		out2 = chooseHere(in2, n1);
% 	else
% 		out1 = chooseHere(in1, n2);
% 		out2 = in2;
% 	end
% end

function out = chooseHere(in, n)
	tmp1 = randperm(length(in));
	tmp2 = sort(tmp1(1 : n));
	out = in(tmp2);
end