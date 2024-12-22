function ME = fcv03_save4parallel( fname, v )
	% % % 13/05/14 written by wp - save within parallel computation
	
	ME = [];
	try
		save(fname, 'v', '-v7.3');
	catch ME
		fprintf(ME.message);
	end

end % end of function

