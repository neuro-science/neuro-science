% % % Filt the data and then rescale
% % % Input and output data are both continuous
% % % Need FiltFiltM package files (+FiltX)
function [data, out] = cmp10_rawmeg_filter(data, srate, freqHighPass, freqLowPass, forder, infoTag, filtOnDim)

% % % rewritten 13/12/2015 by wp and add dimension choice and filter type choice
% % % rewritten 23/07/2014 by wp and remove resampling
% % % updated 5/9/2013 to make it more general 
% % % updated 9/4/2013 to add time estimation

	%% initialization and set parameters
	tic;
	if nargin < 7 || isempty(filtOnDim)
		filtOnDim = 1;
	end
	if nargin < 6 || isempty(infoTag)
		isPlot = 0;
		fid = 1;
	else
		isPlot = 1;
		fid = fopen([infoTag, '.txt'], 'a');
	end
	if nargin < 5
		error('At leaset 5 inputs needed [data, srate, freqHighPass, freqLowPass, forder]!\n');
	end
	
	%% construct and apply filters
	% % % prepare a and b for filter
	if isempty(freqLowPass)
		fprintf(fid, 'HighPass Filter will be used!\n');
		[out.filterB, out.filterA] = butter(forder, 2 * freqHighPass / srate, 'high'); 
	elseif isempty(freqHighPass)
		fprintf(fid, 'LowPass Filter will be used!\n');
		[out.filterB, out.filterA] = butter(forder, 2 * freqLowPass / srate, 'low'); 
	elseif freqHighPass < freqLowPass
		fprintf(fid, 'BandPass Filter will be used!\n');
		[out.filterB, out.filterA] = butter(forder, 2 * [freqHighPass freqLowPass] / srate); 
	else
		fprintf(fid, 'BandStop Filter will be used!\n');
		[out.filterB, out.filterA] = butter(forder, 2 * [freqLowPass freqHighPass] / srate, 'stop'); 
	end
	% % % 3. do filtering
	data = FiltFiltM(out.filterB, out.filterA, data, filtOnDim);	%1st filter for designated range
	
	% % % 4. plot if needed
	if isPlot
		fh = figure('visible','off');
		plot(data);
		print(fh, '-dpng', [infoTag, '_allChans_after_Filter']);
		close(fh);
	end

	%% return data
	fprintf(fid, 'After %4.2f minutes, filtering was done!\n', toc/60);
	if isPlot
		fclose(fid);
	end
end %% end of function