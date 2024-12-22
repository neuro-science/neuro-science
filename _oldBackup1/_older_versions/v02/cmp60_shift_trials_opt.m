function [isSuccess, trIDMatrix, secsCost] = cmp60_shift_trials_opt (numChans, numTrials, secsLimit)
	if nargin < 3
		secsLimit = 60;
	end
	isSuccess = 0;
	trIDMatrix = zeros(numTrials, numChans);
	trIDMatrix(:, 1) = 1 : numTrials;
	tic;
	counter = 1;
	while toc < secsLimit
		tmp2 = 0;
		while ~all(tmp2(:)) && (toc < secsLimit);
			tmp = randperm(numTrials);
			tmp1 = repmat(tmp, [numChans, 1])';
			tmp2 = tmp1 - trIDMatrix;
		end
		counter = counter + 1;
		toc;
		trIDMatrix(:, counter) = tmp;
		if counter > numChans
			isSuccess = 1;
			break;
		end
	end
	secsCost = toc;
end
