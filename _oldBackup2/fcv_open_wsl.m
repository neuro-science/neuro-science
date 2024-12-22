function fcv_open_wsl(fileName, filePath)

	WSL_PATH = '\\WSL$\Ubuntu-22.04\home\pwang\';
	DEFAULT_FOLDER = 'git\weg\4TP\';

	if nargin < 2 || isempty(filePath)
		filePath = DEFAULT_FOLDER;
	end
	
	edit (fullfile(WSL_PATH, filePath, fileName));
end