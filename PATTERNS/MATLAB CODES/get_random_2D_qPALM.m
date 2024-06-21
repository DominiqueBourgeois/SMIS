% Random distribution for qPALM

n=1024; % # of lines
m=1024; % # of columns
I=zeros(n,m);

N=2000; % # of molecules

im_out='C:\Users\bourgeoi\Documents\MATLAB\SMIS\SMIS_vsn2.3\PATTERNS\2D\QPALM\rand_qPALM_1024_2.tif';

xy=rand(2,N);

x=1+round(xy(1,:)*(n-1));
y=1+round(xy(2,:)*(m-1));

for k=1:N
    if I(x(k),y(k))==0
        I(x(k),y(k))=k;
    else
        while I(x(k),y(k))~=0
            x(k)=1+round(rand*(n-1));
            y(k)=1+round(rand*(m-1));
        end
        I(x(k),y(k))=k;
    end
end

imagesc(I)
colormap('gray');
axis image

disp('Writing image ...')
imwrite(uint16(I),im_out);
disp('Done')
