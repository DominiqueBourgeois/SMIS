function     a_h_all = create_pattern(im_par, n_dyes, simul_3D)

n_h=im_par.n*im_par.binning;
m_h=im_par.m*im_par.binning;
nz_h=im_par.nz*im_par.binning;
raster_h=im_par.raster/im_par.binning; % raster size of the high-resolution template image
if simul_3D==0
    a_h = double(zeros(n_h,m_h));
elseif simul_3D==1
    if fix(abs(im_par.pattern_depth)/im_par.raster)>=im_par.nz/2
        error('Increase nz to reach this pattern depth or decrease pattern_depth !');
    end
    a_h = double(zeros(n_h,m_h,nz_h));
end
%Figure out the pattern
%how many pixels is a line space
%     line_space_pixels=fix(im_par.line_space/raster_h);
first_line_pos_pixels=ceil(im_par.first_line_pos/raster_h);
line_thick_pixels=fix(im_par.line_thick/raster_h);
line_increment=fix((im_par.line_thick+im_par.line_space)/raster_h);
last_line_pos_pixels=min(n_h,first_line_pos_pixels+(im_par.number_of_lines-1)*line_increment);
if line_thick_pixels==0
    error(['Thickness of lines must exceed raster size divided by binning factor, ie: ',num2str(im_par.raster/im_par.binning), ' [nm]']);
end

for i=1:line_thick_pixels
    if simul_3D==0
        a_h(first_line_pos_pixels+i-1:line_increment:min(m_h,last_line_pos_pixels+line_thick_pixels-1), : )=100;
    elseif simul_3D==1
        a_h(first_line_pos_pixels+i-1:line_increment:min(m_h,last_line_pos_pixels+line_thick_pixels-1), : ,nz_h/2+fix(pattern_depth*im_par.binning/im_par.raster))=100;
    end
end

size_a_h=size(a_h);
if numel(size_a_h)==3 && simul_3D==0
    disp('Error: The pattern is a 3D pattern ! Only valid for simulations in 3D');
    return;
end
if numel(size_a_h)==2 && simul_3D==1
    disp('Error: The pattern looks like a 2D pattern, but simulation is set in 3D');
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
        if numel(size_a_h)==3; a_h_all(i,:,:,:)=a_h; end
        if numel(size_a_h)==2; a_h_all(i,:,:)=a_h; end
    end
end
clear('a_h');


end