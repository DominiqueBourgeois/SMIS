function     [a_h_all, n,m,nz, n_sp] = read_patterns(in_images_dir, in_images_names, im_par,simul_3D)

% PURPOSE:
% Reads in the patterns
%
% INPUTS:
%	in_images_dir: array of the directory names where the patterns for each dye are
%	in_images_names: array of the pattern names (.tif format for 2D, .mat format for 3D) for each dye are
%	im_par: imaging parameters
%   simul_3D: flag for 3D imaging
%
% OUTPUTS:
%	n,m,nz: size of the patterns. They must all have the same size
%   a_h_all = the patterns for all dyesh
%   n_sp: # of subpatterns for each dye
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, November 2019: reads in # of subpatterns for each dye 

n_dyes=size(in_images_names,1);
n_sp=ones(1, n_dyes); % Array containing the # of subpatterns in each image
[a_h,n,m,nz]=get_pattern(fullfile(char(in_images_dir(1)),char(in_images_names(1))), im_par);
% n_sp(1)=max(a_h(:))-min(a_h(:))+1;
n_sp(1)=numel(unique(a_h));
size_a_h=size(a_h);
if numel(size_a_h)==3 && simul_3D==0
    disp('Error: The pattern is a 3D pattern ! Only valid for simulations in 3D: check "simul_3D or change pattern files !');
    a_h_all=[];
    return;
end
if numel(size_a_h)==2 && simul_3D==1
    disp('Error: The pattern looks like a 2D pattern, but simulation is set in 3D:  check "simul_3D or change pattern files !');
    a_h_all=[];
    return;
end
if numel(size_a_h)==3 && simul_3D==1
    a_h_all=zeros(n_dyes,size_a_h(1),size_a_h(2),size_a_h(3));
    a_h_all(1,:,:,:)=a_h;
end
if numel(size_a_h)==2 && simul_3D==0
    a_h_all=zeros(n_dyes,size_a_h(1),size_a_h(2));
    a_h_all(1,:,:)=a_h;
end

if n_dyes>1
    for i=2:n_dyes
        [a_h,n,m,nz]=get_pattern(fullfile(char(in_images_dir(i)),char(in_images_names(i))), im_par);
        if size_a_h~=size(a_h)
            disp(['Input pattern number: ',num2str(i),' has a different size than pattern number 1']);
            return
        end
        n_sp(i)=numel(unique(a_h));
        if numel(size_a_h)==3; a_h_all(i,:,:,:)=a_h; end
        if numel(size_a_h)==2; a_h_all(i,:,:)=a_h; end
    end
end

for k=1:n_dyes
    disp(['# of subpatterns for dye number: ', num2str(k),  ' : ', num2str(n_sp(k))]);
end
clear('a_h');

end