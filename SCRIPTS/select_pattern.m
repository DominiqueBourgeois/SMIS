function     a=select_pattern(a_all,i)

s=size(a_all);
d=numel(s)-1;
if d==1 % 1D matrix
    a=squeeze(a_all(i,:));
elseif d==2 % 2D matrix
    a=squeeze(a_all(i,:,:));
elseif d==3 % 3D matrix
    a=squeeze(a_all(i,:,:,:));
end
end