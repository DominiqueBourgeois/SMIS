function sm_par = get_velocity_radii(im_par, sm_par)

% PURPOSE:
% Get velocity radii for sm's diffusing in 2D or 3D at
% speed V
%
% INPUTS:
%	im_par: the imaging parameters
%	sm_par: the sm parameters
%
% OUTPUTS:
%	sm_par: the updated sm parameters
%
% MODIFICATION HISTORY:
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3
%	D.Bourgeois, January 2021: Corrected bug

r1=zeros(numel(sm_par.V),1);
r2=zeros(numel(sm_par.V),1);

for i=1:numel(sm_par.V)
    if sm_par.V(i)>0    
        %Define radius for actual distance covered at speed V in
        %frametime + addtime
        r1(i)=sm_par.V(i)*(im_par.frametime+im_par.addtime)/im_par.raster*im_par.binning;
        
        %Define radius for distance corresponding to persistence length
        %for determination of initial velocity directions
        if sm_par.persistence_length(i)==-1 % Handle automatic determination of persistence length based on chosen velocities
            sm_par.persistence_length(i)=1e-3*sm_par.V(i)*(im_par.frametime+im_par.addtime);
            disp(['Persistence length for velocity ',num2str(sm_par.V(i)),' set to: ',num2str(sm_par.persistence_length(i)),' um']);
        end        
        r2(i)=1e+3*sm_par.persistence_length(i)/im_par.raster*im_par.binning;
    end
end


%2D case
if im_par.simul_3D==0
    for i=1:numel(sm_par.V)
        if sm_par.V(i)>0
            [sm_par.V_circle(i).x,sm_par.V_circle(i).y]=get_circle_coordinates(r1(i)); % Specific to velocity of molecules
            [sm_par.V_init_dir(i).x,sm_par.V_init_dir(i).y]=get_circle_coordinates(r2(i)); % Specific to persistence length 
        end
    end
elseif im_par.simul_3D==1 
    for i=1:numel(sm_par.V)
        if sm_par.V(i)>0
            [sm_par.V_circle(i).x,sm_par.V_circle(i).y,sm_par.V_circle(i).z]=get_sphere_coordinates(r1(i)); % Specific to velocity of molecules
            [sm_par.V_init_dir(i).x,sm_par.V_init_dir(i).y,sm_par.V_init_dir(i).z]=get_sphere_coordinates(r2(i)); % Specific to persistence length 
        end
    end
end
