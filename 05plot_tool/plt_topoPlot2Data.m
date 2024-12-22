function [pZ, pX, pY] = plt_topoPlot2Data(Z, posiChan, XYRange, nPixels1D)
% % % updated 25/08/2014 - data only
% % % updated 16/05/2013 - improve figure property
% % % updated 15/05/2013 - add central sensor bold
% % % topo plot done ?

	%% 1. para set
	% % % set some parameters
	if nargin < 4
		nPixels1D = 1000;
	end
	if nargin < 3
		XYRange = 0.5;
	end
	
	% % % prepare data xy
	tmp = linspace(-XYRange, XYRange, nPixels1D);
	[pX, pY] = meshgrid(tmp, tmp);

	% % % prepare data z
	if verLessThan('matlab', '8.1')
		F = TriScatteredInterp(posiChan, Z);
	else
		F = scatteredInterpolant(posiChan, Z);
	end
	pZ = F(pX, pY);

	% % % Remove outsiders
	outSiders = (sqrt(pX.^2 + pY.^2) > XYRange);
	pZ(outSiders) = NaN;

end %end of function