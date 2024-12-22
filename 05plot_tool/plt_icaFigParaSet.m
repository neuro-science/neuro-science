function fp = plt_icaFigParaSet(varargin)

	%% 1. para for figure and axis positions etc
	% % % figure properties	
	fp.fColor = [1 1 1];
	fp.fOri = 'Portrait';
% 	fp.fOri = 'Landscape';
% % % 	% % % figure size, about A4 ratio in 1080p screen	
% % % 	fp.fStartXPx = 50;
% % % 	fp.fStartYPx = 50;
% % % 	fp.fWidthPx = 1200;
% % % 	fp.fHeightPx = 850;

	% % % 	figure size, keep A4 ratio for printing
	fp.fgSz = round([0 0 11.7 8.27] * 300);
	
	% % % border settings	
	fp.fStartX = 0.03;	%whole figure start point in X
	fp.fEndX = 0.97;
	fp.fStartY = 0.03;
	fp.fEndY = 0.97;
	fp.bMidX1 = 0.35;	%border of left and right part in X
	fp.bMidX2 = 0.55;	%border of left and right part in X
	fp.bTopY2 = 0.85;
	fp.bTopY = 0.9;
	fp.bBotY = 0.55;
	fp.spAX1D = 0.03;
	
	%% 2. para for text
	fp.txColorAX	= [1 1 1];
	fp.txRangeX = [0, 1];
	fp.txRangeY = [0, 1];
	fp.txStartX = 0;
	fp.txStartY = 0.5;
	fp.txS = '\n[Here is the Position for component name tag]';
	fp.txS2 = '[Here is the Position for correlations]';	
	fp.txSZ = 12;
	fp.txFT = 'Helvetica';
	
	%% 3. curve/trl plot para
	fp.cvColorAX	= [1 1 1];
	fp.cvNum = 10;
	fp.cvRangeX = [0, 1];
	fp.cvRangeY = [0, 10];
	fp.cvTmpPtsNum = 2401;
	fp.cvColor = [0 0 1];
	fp.cvWidth = 0.1;
	fp.trColorAX = [1 1 1];
	
	%% 4. topo plot
	fp.tpColorAX	= [1 1 1];
	fp.tpTmpPtsNum = 1000;
	fp.tpRangeXY = 0.5;
	fp.tpAspRatio = [255 255 1];
	
	%% 5. power plot
	fp.pwColorAX	= [1 1 1];
	fp.pwTmpPtsNum = 100;
	fp.pwRangeX = [0, 200];
	fp.pwRangeY = [0, 1];
	fp.pwColor = [1 0 0];
% 	fp.pwWidth = 2;
	
	%% 6. use iput para to replace default para
	if (rem(length(varargin), 2))
		error('Optional parameters should always go by pairs');
	else
		fs = fieldnames(fp);
		for k = 1 : 2 : (length(varargin)-1)
			if ~ischar (varargin{k}),
				error (['Unknown type of optional parameter name (parameter' ...
				' names must be strings).']);
			else
				yesFieldID = strcmpi(varargin{k}, fs);
				if any(yesFieldID)
					cmd = ['fp.', fs{yesFieldID}, '= varargin{k + 1};'];
					eval(cmd);
				else
					error ('Unknown field names, check again!');
				end
			end
		% change the value of parameter
		end
	end
	
	%% 7. induced parameters
	% % % axes positions
	fp.spAXPos = [fp.spAX1D, fp.spAX1D, -fp.spAX1D, -fp.spAX1D];
	fp.txPos	= [fp.fStartX, fp.bTopY, fp.fEndX - fp.fStartX, fp.fEndY - fp.bTopY] + fp.spAXPos; %text axes position
	fp.txPos2	= [fp.fStartX, fp.bTopY2, fp.fEndX - fp.fStartX, fp.bTopY - fp.bTopY2] + fp.spAXPos; %text axes position
	fp.cvPos	= [fp.fStartX, fp.fStartY, fp.fEndX - fp.fStartX, fp.bBotY - fp.fStartY] + fp.spAXPos; %curve plot axes position
	fp.tpPos	= [fp.fStartX, fp.bBotY, fp.bMidX1 - fp.fStartX, fp.bTopY2 - fp.bBotY] + fp.spAXPos; %topo plot axes position
	fp.pwPos	= [fp.bMidX2, fp.bBotY, fp.fEndX - fp.bMidX2, fp.bTopY2 - fp.bBotY] + fp.spAXPos; %power plot axes position
	fp.trPos	= [fp.bMidX1, fp.bBotY, fp.bMidX2 - fp.bMidX1, fp.bTopY2 - fp.bBotY] + fp.spAXPos; %trl plot axes position
	
	% % % curve data induction
	fp.cvDataX = linspace(fp.cvRangeX(1), fp.cvRangeX(2), fp.cvTmpPtsNum)';
	tmp = linspace(fp.cvRangeY(1), fp.cvRangeY(2), fp.cvNum + 1);
	tmp(end) = [];
	fp.cvDataY = bsxfun(@plus, (1 + sin(bsxfun(@plus, fp.cvDataX, rand(1, fp.cvNum)) * 20 * pi)) / (2 * fp.cvNum), tmp);
	fp.cvXTick = fp.cvDataX([1, end]);
	fp.cvXTickLabel = num2str(fp.cvXTick);
	fp.cvYTick = linspace(fp.cvRangeY(1), fp.cvRangeY(end), fp.cvNum + 1);
% 	fp.cvYTickLabel = num2str(fp.cvYTick');
	fp.cvYTickLabel = [];

	% % % topo data induction
	fp.tpXY = linspace(-fp.tpRangeXY, fp.tpRangeXY, fp.tpTmpPtsNum);
	fp.tpRangeX = [-fp.tpRangeXY, fp.tpRangeXY];
	fp.tpRangeY = [-fp.tpRangeXY, fp.tpRangeXY];
	% % % 	[fp.tpX, fp.tpY] = meshgrid(tmp, tmp);
	% % % 	fp.tpZ = 1 ./ (sqrt((fp.tpX - fp.tpTmpPtsNum/2).^2 + (fp.tpY - fp.tpTmpPtsNum/2).^2)/fp.tpTmpPtsNum + 1) - 0.5;
% 	fp.tpZ = (fp.tpX - 0.5).^2 + (fp.tpY - 0.5).^2;
	fp.tpZ = bsxfun(@plus, fp.tpXY.^2, (fp.tpXY').^2);

	% % % power data induction
	fp.pwX = linspace(fp.pwRangeX(1), fp.pwRangeX(2), 100)';
	fp.pwY = rand(fp.pwTmpPtsNum, 1) / 3 + (sin(fp.pwX * 4 * pi) + 1) / 3;
	
	% % % trl data
	NN = 1000;
	fp.trx = rand(NN, 1);
	fp.try = 1 : NN;
	
end % end of function

