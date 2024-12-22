function [pac, d] = cmp53_cross_freq (cfg, data, freqs)
% % % 22/01/16	written by wp: 
% % % 	data: nTimePoints x ... x ...
% % %		freqs: filter ranges

	%% 1. check inputs
	% % % data format	
	data = double(permute(data, [2 1 3])); %expect as [ch, pnt, tr] -> [pnt ch tr]
	[nps, nchs, ntrs] = size(data);
	% % % config paras
	if ~isfield(cfg, 'fOrder')
		cfg.fOrder = 4;
	end
	if ~isfield(cfg, 'sRate')
		cfg.sRate = 1000;
	end
	if ~isfield(cfg, 'tID')
		cfg.tID = 1:nps;
	end
	
	%% 2. prepare data
	% % % prepare filter
	nfs = size(freqs, 1);
	if size(freqs, 2) ~=4
		error('The frequency shall be [lowHighpass, lowLowpass, highLowpass, highLowpass] format!');
	end
	% % % get data
	d = zeros([nps, nchs, ntrs, nfs, 2]); %[pnt ch tr fpair l/h]
	for k = 1 : nfs
		[fB{k, 1}, fA{k, 1}] = butter(cfg.fOrder, 2 * [freqs(k, 1) freqs(k, 2)] / cfg.sRate); 
		[fB{k, 2}, fA{k, 2}] = butter(cfg.fOrder, 2 * [freqs(k, 3) freqs(k, 4)] / cfg.sRate); 
		d(:, :, :, k, 1) = cmp08_preCrossFrequency(data, fA{k, 1}, fB{k, 1}, 'l');
		d(:, :, :, k, 2) = cmp08_preCrossFrequency(data, fA{k, 2}, fB{k, 2}, 'h');
	end
	
	
	%% 3. work on coupling
	pac = zeros(nchs, nchs, ntrs, nfs);
	for k1 = 1 : nfs
		for k2 = 1 : ntrs
			pac(:, :, k2, k1) = corr(d(cfg.tID, :, k2, k1, 1), d(cfg.tID, :, k2, k1, 2));
		end
	end
		
end %end of function

