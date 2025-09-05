function MySample=get_SMIS_2DSquaresCircles(par)


n=par.x_dim; %number of column
m=par.y_dim; %number of slices

is_square=par.is_square;

qPALM_option=par.qPALM;

%spacing between lines [pixel]
s=par.s; % size or diameter
% N=par.N; % # of elements

n_row=par.n_row;
n_col=par.n_col;

r_shift=par.random_shift;

%safety border
border=par.border;

%Actual values to use
n2=n-border;
m2=m-border;

feature_id=1;
%%
MySample = double(zeros(m,n));

if (n_col>n_row && m<n) || (n_col<n_row && m>n)
    tmp=n_row;
    n_row=n_col;
    n_col=tmp;
end

inc_r=fix(n2/n_row)+1; % # of increments along rows
inc_c=fix(m2/n_col)+1; % # of increments along columns

x0=inc_r:inc_r:inc_r*n_row;
y0=inc_c:inc_c:inc_c*n_col;

offset_x=(x0(1)+x0(end)-(n2-1))/2;
offset_y=(y0(1)+y0(end)-(m2-1))/2;

x0=x0-fix(offset_x)+border;
y0=y0-fix(offset_y)+border;

v=1; % pattern id, to be increased for qPALM

[x,y]=meshgrid(1:n,1:m);
disp('Creating pattern ...')
for i=1:n_row
    x_c=round(x0(i));
    for j=1:n_col
        y_c=round(y0(j));

        if r_shift>0
            x_c=x_c+r_shift*(rand-0.5);
            y_c=y_c+r_shift*(rand-0.5);
        end

        if is_square==1
            if mod(s,2)==1
                t=abs(x-x_c)<=fix(s/2) & abs(y-y_c)<=fix(s/2);
            else
                t=x<(x_c+s/2+1) & x>=(x_c-s/2+1) & y<(y_c+s/2+1) & y>=(y_c-s/2+1);
            end
        else
            t=((x-x_c).^2+(y-y_c).^2)<=s^2;
        end

        if qPALM_option==1
            MySample(t==1)=v;
            v=v+1;
        else
            MySample=MySample | t;
        end
    end
end
disp('Done !')

MySample=MySample*feature_id;

% if border>0
%     MySample(1:border,:)=0;
%     MySample(end-border:end,:)=0;
%     MySample(:,1:border)=0;
%     MySample(:,end-border:end)=0;
% end

% Make sure qPALM patterns start at 1
if qPALM_option==1
    u_val=unique(MySample);
    if (numel(u_val)-1)~=n_row*n_col
        disp(['Number of clusters created: ',num2str(numel(u_val)-1), ' not equal to requested number (',num2str(n_row*n_col),') !']);
        %Reorder the clusters from 1 to numel(u_val)-1
        % Get the unique pixel values and sort them
        disp('Reassigning pixels ...');
        sortedValues = sort(u_val);

        % Create a mapping from the original values to the new values
        valueMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
        for i = 1:length(sortedValues)
            valueMap(sortedValues(i)) = i-1;
        end

        % Reassign the pixel values
        newImage = zeros(size(MySample));
        for i = 1:numel(MySample)
            newImage(i) = valueMap(MySample(i));
        end
        MySample=newImage;
    end
end

%%
disp('Done !');


%% Show the cell
figure(1)
clf
set(gcf,'Color','w')
imagesc(MySample);
axis image
colormap('gray')
xlabel('X [pixel]')
ylabel('Y [pixel]')
if is_square==1
    title('Squares')
else
    title('Circles')
end

