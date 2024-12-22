function [tf_data, myWinSizes, nTapers, myWins, myTapers, myIndFFT] = ...
	mat_multitaperFFT(data, sRate, timeStart, xTimes, yFreqs, winSizes, resFreq, fid)
% % % 16/02/24 updated by wp: winsize option may be vector or scalar
% % % 09/02/24 updated by wp: more possiblities of data input format
% % % 11/07/17 updated by wp: newer filename
% % % 22/11/16 updated by wp: adjust the taper numbers for small window
% % % 01/10/15 updated by wp: myWins initialize correct dimension
% % % 23/07/15 updated by wp: power only output
% % % 25/02/15 updated by wp: remove outputs in some cases
% % % 26/06/14	updated by wp: rename for newer format
% % % 12/06/14 updated by wp: direct output to fid if available

% % % tf_data = wpa_0Gen_F100MultitaperFFT(data, sRate, timeStart, xTimes, yFreqs, winSizes, resFreq)
% % % The data shall be in channel x point x trial format
% % % The unit of timeStart, xTimes and winSize is ms and the unit of yFreq is Hz
% % % frequency resolution shall be one number indicating octaves(log2(fMax/fMin))
% % % winSize shall be one number in ms, later may be extended
% % % output is a cell array tf_data{time, freq}, the element is [channel, trial, taper]

	%% preparison
	tic;
	if nargin < 8 || isempty(fid)
		fid = 1;
		msflag = false;	%do log? - default no
		pflag = false;		%do power only? - default no
		tflag = true;		%data format, default [channel x point x trial] later to [points, trial*channel]
	elseif fid == -1
		msflag = false;
		pflag = true;
		tflag = true;
	elseif fid == -2
		msflag = false;
		pflag = true;
		tflag = false;	%alternative input is [point x trial x channel]
	elseif fid == -3
		msflag = false;
		pflag = false;
		tflag = false;
	else
		msflag = true;
		pflag = false;
		tflag = true;
	end
	
	% % % get data size and reshape data to [points, trial*channel]
	if tflag
		[nChans, nPoints, nTrials] = size(data);
		data = reshape(permute(data, [2 3 1]), [nPoints, nTrials * nChans]);
	else	%alternative input is [point x trial x channel]
		[nPoints, nTrials, nChans] = size(data);
		data = reshape(data, [nPoints, nTrials * nChans]);
	end
	
	% % % data info for tf request
	nTimes = length(xTimes);
	nFreqs = length(yFreqs);
	
	% % % output definition
	if pflag
		tf_data = zeros(nTrials, nChans, nTimes, nFreqs);
	else
		tf_data = cell(nTimes, nFreqs);
	end
	
	% % % calculate time indices to check data boundary
	xPoints = round((xTimes - timeStart) * sRate /1000) + 1;
	if min(xPoints) < 1
		fprintf('%d is requested, ', min(xPoints));
		error('time range was exceeded in the left!');
	end
	if max(xPoints) > nPoints
		fprintf('%d vs %d, ', max(xPoints), nPoints);
		error('Time range was exceeded in the right!');
	end
	
	% % % the window sizes will change slightly to meet integer number of tapers
	if numel(winSizes) == 1 %uniform size extend to vector
		winSizes = zeros(nFreqs, 1) + winSizes;
	elseif numel(winSizes) ~= nFreqs
		error('Size error of input window size!');
	end
	myWinSizes = winSizes / 1000 * sRate;
	myWins = zeros(nFreqs, 2);
	nTapers = zeros(nFreqs, 1);
	myTapers = cell(nFreqs, 1);
	myIndFFT = zeros(nFreqs, 1);
	
	%% taper calculation
	for f = 1 : nFreqs
		df = yFreqs(f) * (2.^(resFreq/2) - 2.^(-resFreq/2));
		if yFreqs(f) >= 15.9 && df * winSizes(f) > 2000	% 16+
% 		theFreq = yFreqs(f);
			nTapers(f) = round(df * winSizes(f) /1000) - 1;
		else
			nTapers(f) = 1;
		end
		myWinSizes(f) = round((nTapers(f) + 1) / df * sRate);
		myWins(f, 1) = -floor(myWinSizes(f)/2);
		myWins(f, 2) = myWins(f, 1) + myWinSizes(f) - 1;
		tmpFreq = linspace(0, sRate/2, ceil(myWinSizes(f)/2));
		[tmp, myIndFFT(f)] = min(abs(tmpFreq - yFreqs(f)));
		[e, v] = dpss(myWinSizes(f), myWinSizes(f) * df / sRate / 2);
		myTapers{f} = permute(e(:, 1 : nTapers(f)), [1 3 2]);	%facilate bsxfun
		clear tmp tmpFreq df e v;
	end

	%% calculation
	for f = 1 : nFreqs
		for t = 1 : nTimes
			% % % select data
			sid = xPoints(t) + (myWins(f, 1) : myWins(f, 2));	%indices needed
			padFlag = 0;	%initialize
			n1 = 0;
			n2 = 0;
			id1 = find(sid < 1);	%left exceed
			if ~isempty(id1)
				n1 = length(id1);
				if msflag
					fprintf(fid, 'Indices exceeded left side in [%d ms,%6.2f Hz], padding %d zeros...\n', ...
						xTimes(t), yFreqs(f), n1);
				end
				sid(id1) = 1;	%use first points
				padFlag = 1;
			end
			id2 = find(sid > nPoints);	%right exceed
			if ~isempty(id2)
				n2 = length(id2);
				if msflag
					fprintf(fid, 'Indices exceeded right side in [%d ms,%6.2f Hz], padding %d zeros...\n', ...
						xTimes(t), yFreqs(f), n2);
				end
				sid(id2) = nPoints;	%use last points
				padFlag = 1;
			end
			theData = data(sid, :);	%get data
			if padFlag
				theData([1 : n1, end - n2 + 1 : end], :) = 0; %pad zeros
			end
			% % % apply taper [points, trial*channel, taper]
			theData = bsxfun(@times, theData, myTapers{f});
			% % % do FFT			
			theTF = fft(theData);
			if pflag
				tmp = nanmean(theTF(myIndFFT(f), :, :) .* conj(theTF(myIndFFT(f), :, :)), 3);
				tf_data(:, :, t, f) = reshape(tmp, [nTrials, nChans]);
			else
				% % % sort data to [channel, trial, taper]			
				tf_data{t, f} = permute(reshape(theTF(myIndFFT(f), :, :), ...
					[nTrials, nChans, nTapers(f)]), [2 1 3]);
			end
% 			clear theTF theData padFlag sid id1 id2 n1 n2;
			fprintf('.');
		end
		fprintf('o\n');
		if msflag
			fprintf(fid, 'Done for frequency %6.2f after %6.1f seconds. \n', yFreqs(f), toc);
		end
	end
end %end of function