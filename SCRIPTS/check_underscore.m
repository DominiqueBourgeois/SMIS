function my_new_str = check_underscore(my_str)
% Add a backlash to string if it has a underscore character !
    my_new_str = strrep(my_str,'_','\_');
end