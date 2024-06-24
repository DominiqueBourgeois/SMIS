function out=on2one(in,str1,str2)

if isnumeric(in)
    if in==0
        out=str1;
    else
        out=str2;
    end
else
    if strcmp(in,str1)==1
        out=0;
    else
        out=1;
    end
end

