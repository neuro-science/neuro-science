function out = plt_stringConverter4ProperDisplay(in)
% % % written by wp 04/04/2023 to convert strings 
% % % such that Tex interpreter would not eat the characters
   out = strrep(in, '\', '\\');
   out = strrep(out, '_', '\_');
   out = strrep(out, '^', '\^');
end