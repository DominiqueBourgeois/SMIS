function r = check_dir(my_dir)
% Check if a directory exists !
if isequal(exist(my_dir, 'dir'),7) % 7 = directory.
    % We have a folder!
    r=1;
else
    % We have an invalid file or folder.
    r=0;
    warning('Directory does not exist!');
end