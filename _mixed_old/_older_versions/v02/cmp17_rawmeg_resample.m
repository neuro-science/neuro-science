function data = cmp17_rawmeg_resample(data, ratio)
% % % 28/07/14	written by wp: no filter included !!!
	data = data(1 : ratio : end);

end %end of function

