function s = fcv05_rename_fields (s, str1, str2)
% % % rename fields of a structure: str1 to str2
% % % 23/05/2014 by wp
	sf = fieldnames(s);
	for k = 1 : numel(sf)
		eval(['tmp = s.', sf{k}, ';']);
		if ischar(tmp) 
			tmp = strrep(tmp, str1, str2);
			fprintf('   replaced in field: %s\n', sf{k});
		elseif iscell(tmp)
			try
				tmp = strrep(tmp, str1, str2);
				fprintf('   replaced in field: %s\n', sf{k});
			catch ME
				fprintf('cell array but not strings: %s\n', sf{k});
			end
		elseif isstruct(tmp)
			tmp = fcv05_rename_fields (tmp, str1, str2);
			fprintf('Go deeper in field: %s\n', sf{k});
		end
		eval(['s.', sf{k}, ' = tmp;']);
	end %end of fields loop
			

end