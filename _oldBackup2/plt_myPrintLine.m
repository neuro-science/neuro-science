function plt_myPrintLine(str, varargin)
% % % [*1] modified 06/10/2022 such that the default one return not counted
% % % First written early 2022 for the MEG experiment


	% % % security check
	if nargin < 1
		error('You need to provide the string to print!');
	end
	% % % default parameters
	label = '-';	%The charactor serve as lables
	% % % [*1] below modified 06/10/2022 such that the default one return not counted
	N = [0 0];		%number of returns before and after this line
	L = 36;			%Half of default line width
	
	% % % evaluate the parameters from inputs	
	idx = find(strcmp('label', varargin));
	if ~isempty(idx)
		label = varargin{idx + 1};
		varargin([idx, idx + 1]) = [];
	end
	
	idx = find(strcmp('L', varargin));
	if ~isempty(idx)
		L = varargin{idx + 1};
		varargin([idx, idx + 1]) = [];
	end
	
	idx = find(strcmp('N', varargin));
	if ~isempty(idx)
		N = varargin{idx + 1};
		varargin([idx, idx + 1]) = [];
	end
	
	% % % do the job
	numS = floor(L - numel(str)/2);
	for k = 1 : N(1)
		fprintf('\n');
	end
	fprintf([repelem(label, numS), str, repelem(label, numS)], varargin{:});
	% % % [*1] below modified 06/10/2022 such that the default one return not counted
	for k = 1 : N(2) + 1
		fprintf('\n');
	end
end