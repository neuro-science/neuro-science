% % % written by wp 11/12/2023 to integrate previous plot functions for 2-D topo plot

function [scaledDataRGB, x, y] = plt_topoEEG(Z, loc2d, axHandle, zThreshold, maxValue4Z, ...
	xRange, yRange, pointsPerUnit, colorStyle)
% % % plt_topoEEG: Plots a 2D topographic map of EEG data.
% % %
% % % Inputs:
% % %   Z					- EEG data at each sensor (num_of_sensor x 1).
% % %   loc2d				- Sensor locations (num_of_sensor x 2) indicating x, y coordinates.
% % %   axHandle		- handle for the figure.
% % %   pointsPerUnit	- Resolution (number of pixels per unit of range).
% % %   colorStyle		- Style of the colormap (e.g., 'jet', 'hot', 'custom').
% % %   xRange				- Range for x-axis, [minX maxX]. If empty, computed from sensor.
% % %   yRange				- Range for y-axis, [minY maxY]. If empty, computed from sensor.
% % %   zThreshold		- Threshold for color scaling.
% % %   maxValue4Z		- Maximum absolute value for color scaling.
% % % Outputs:
% % %   scaledDataRGB	- plotting data at each point of the figure.
% % %   x					- x data
% % %   y					- y data

    % Set default values and compute ranges if not provided
	if nargin < 9, colorStyle = 'jet'; end
	if nargin < 8, pointsPerUnit = 1001; end
	if nargin < 7 || isempty(yRange), yRange = [min(loc2d(:,2)) max(loc2d(:,2))]; end
	if nargin < 6 || isempty(xRange), xRange = [min(loc2d(:,1)) max(loc2d(:,1))]; end
	if nargin < 5, maxValue4Z = max(abs(Z)); end
	if nargin < 4, zThreshold = 0; end
	if nargin < 3 || isempty(axHandle), axHandle = false; end
	if nargin < 2, error('We need at least the data and location of sensors!'); end

    % Prepare data grid
    x = linspace(xRange(1), xRange(2), pointsPerUnit);
    y = linspace(yRange(1), yRange(2), pointsPerUnit);
    [pX, pY] = meshgrid(x, y);

    % Data interpolation
    F = scatteredInterpolant(loc2d, Z, 'natural', 'none');
    pZ = F(pX, pY);

    % Ellipse parameters (modify these as needed)
    ellipseCenter = [(xRange(1) + xRange(2)) / 2, (yRange(1) + yRange(2)) / 2];
    ellipseAxes = [diff(xRange) / 2, diff(yRange) / 2]; % Semi-major and semi-minor axes

    % Set outsiders to NaN based on ellipse
    outSiders = ((pX - ellipseCenter(1)).^2 / ellipseAxes(1)^2 + ...
		 (pY - ellipseCenter(2)).^2 / ellipseAxes(2)^2) > 1;
    pZ(outSiders) = NaN;

    % Scale and color data
    scaledData = scaleAndColorData(pZ, zThreshold, maxValue4Z, colorStyle);
    
    % Set NaNs to white
    nanMask = isnan(scaledData);
    scaledDataRGB = ind2rgb(gray2ind(scaledData, 256), colormap(colorStyle));
    white = [1 1 1];
    scaledDataRGB(repmat(nanMask, [1 1 3])) = repmat(white, sum(nanMask(:)), 1);

    % Plotting
	 if axHandle
	    im = image(axHandle, x, y, scaledDataRGB);
	    axis equal tight;
	    colormap(colorStyle);
	    colorbar;
	 end
	 
	 % clean up
	 clear white nanMask scaledData outSiders ellipseAxes ellipseAxes F pX pY pZ;
end

function cout = scaleAndColorData(data, zThreshold, maxValue4Z, colorStyle)
    % Scale data and apply color mapping
    data = data / maxValue4Z;
    data(abs(data) > 1) = sign(data(abs(data) > 1));
    data(abs(data) < zThreshold) = 0;

    % Apply custom color scaling if needed
    if strcmpi(colorStyle, 'custom')
        cout = customColorMapping(data); % Custom color mapping function
    else
        cout = data;
    end
end

function cout = customColorMapping(data)
    % customColorMapping: Example of a custom color mapping function.
    % This function can be modified to create a custom color scale.
    %
    % Inputs:
    %   data - Scaled EEG data
    %
    % Output:
    %   cout - Color-mapped data

    % Example of a simple custom mapping: linear gray scale
    cout = (data + 1) / 2; % Map data to [0, 1]
    cout = repmat(cout, [1, 1, 3]); % Replicate for RGB
end

