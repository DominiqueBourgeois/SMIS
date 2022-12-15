function ok = check_input(dir, file)
% PURPOSE:
%   Check if a file exists
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2012.
%

ok=1;
my_file=fullfile(dir,file);
if ~exist(my_file,'file')
    disp(['File ',my_file,' does not exist !']);
    ok=0;
end
end