function [x,y,z,a_sm,ok] = place_single_multimer_molecules_simple(a,im_par,sm_par,n_clusters)
% NAME:
%	PLACE_SINGLE_MOLECULES
%
% PURPOSE:
%	Get xyz coordinates of SMs following a pattern on an image
%
% CATEGORY:
%	Single molecule.
%
% CALLING SEQUENCE:
%	[x,y,z, a_sm] = place_single_molecules(a,par,channel)
%
% INPUTS:
%	a: the image onto place the SMs
%	im_par: some general parameters
%   sm_par: single molecule parameters
%
% OUTPUTS:
%	x = the x coordinates [in raster units]
%   y = the y coordinates [in raster units]
%   z = the z coordinates [in raster units]
%   a_sm = the image, where the SMs have been positionned
%   ok = 1 if execution was okay, 0 if an error occurred
% COMMON BLOCKS:
%	None.
%
% SIDE EFFECTS:
%	None.
%
% RESTRICTIONS:
%	None.
%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011.
%	D.Bourgeois, May 2015: added z coordinate for 3D.
%	D.Bourgeois, March 2021: Added error handling

ok=1;
max_a=double(max(a(:))); % Maximum amplitude in a
if max_a ~= n_clusters
    disp(['Fluorophore: ',sm_par.fluorophore_name]);
    errordlg(['The pattern contains a maximum of: ', num2str(max_a), ' clusters and is not compatible with distribution of dyes within clusters !'],'SMIS Error');
    uiwait
    ok=0;
    return
end

x=[];
y=[];
z=[];
a_sm=a;
a_sm(:)=0; % reset

n_mol_p=sm_par.n_mol*sm_par.n_sp_fr; % # of molecule to place on pattern
n_mol_op=sm_par.n_mol*(1-sm_par.n_sp_fr); % # of molecule to place out of pattern


n_mol_per_cluster=fix(n_mol_p/n_clusters);

if im_par.simul_3D==0 % 2D mode
    for i=1:max_a
        if max_a>1 && fix(i/100) ~=0 && (i/100==fix(i/100))
            clc;
            disp(['Placed molecule [%]: ',num2str(100*i/max_a)]);
        end
        %indices of pixels that might be labeled
        w_a=find(a==i);
        
        %each pixel should be subdivided to account for nearest neighbour distance
        subdiv_factor=max([1,fix(im_par.raster/sm_par.near_dist/im_par.binning)]);
        sub_a=zeros(subdiv_factor,subdiv_factor); % create sub image
        %so that # of subpixels that can potentially be labeled is:
        n_to_label=size(w_a,1)*subdiv_factor^2;
        
        if n_to_label<n_mol_per_cluster
            disp(['Fluorophore: ',sm_par.fluorophore_name]);
            disp(['Number of possibly labeled pixels: ',num2str(n_to_label)]);
            errordlg('Number of molecules to label is too high; retry or reduce # of molecules','SMIS Error');
            uiwait
            ok=0;
            return
        end
        
        %pick randomly locations where labeling molecules will end up
        r=rand(1,n_to_label); %uniformly distributed random number
        [~, sorted_i]=sort(r,'ascend'); %sorted_i are the indices of sorted_r
        sorted_i=sorted_i'; %make it a line vector
        % below: corrected bug 22/04/2021
        sorted_i_a = fix((sorted_i-1)/subdiv_factor^2); %the index of the big pixel in  a 
        sorted_i_r = fix(sorted_i-sorted_i_a*subdiv_factor^2); %the index of the small pixel in sub_a
        sorted_i_a=sorted_i_a(1:n_mol_per_cluster)+1;
        sorted_i_r=sorted_i_r(1:n_mol_per_cluster);
        if max(sorted_i_a) > length(w_a)         
            SMISMessage1=['Fluorophore: ',sm_par.fluorophore_name];
            SMISMessage2='Decrease min distance between fluorophores !';
            SMISMessage=[SMISMessage1 newline SMISMessage2];
            disp(SMISMessage);
            warndlg(SMISMessage,'Warning')
            uiwait
            ok=0;
            return
        end
        
        i_selected=w_a(sorted_i_a);
        for j=1:n_mol_per_cluster
            a_sm(i_selected(j))=a_sm(i_selected(j))+100; %Loop to figure out if multiple SM in 1 pixel
        end
        [x0,y0]=wheresub(i_selected,a); % x is the row #; y is the column #
        [sub_x,sub_y]=wheresub(sorted_i_r,sub_a);
        % Get the xy position of all molecules and add some random components
        new_x=x0-0.5+1/subdiv_factor*((sub_x-0.5)+(rand(n_mol_per_cluster,1)-0.5)); % molecule in pixel x might take positions from x-0.5 to x+0.5
        new_y=y0-0.5+1/subdiv_factor*((sub_y-0.5)+(rand(n_mol_per_cluster,1)-0.5)); % molecule in pixel x might take positions from y-0.5 to y+0.5
        x=[x; new_x];
        y=[y; new_y];
    end
    
    % Place the rest of the molecules out of the pattern
    %indices of pixels that might be labeled
    if n_mol_op > 0
        w_a=find(a==0);
        [x_op,y_op,~,a_sm,ok] = add_single_molecules_on_pattern(a_sm,w_a, sm_par, im_par,n_mol_op);
        if ~ok; return; end
        x=[x; x_op];
        y=[y; y_op];
    end
    z=0*x+1; % Set z to 0
    
elseif im_par.simul_3D==1 % 3D mode
    for i=1:max_a
        %indices of pixels that might be labeled
        w_a=find(a==i);
        
        %each pixel should be subdivided to account for nearest neighbour distance
        subdiv_factor=max([1,fix(im_par.raster/sm_par.near_dist/im_par.binning)]);
        sub_a=zeros(subdiv_factor,subdiv_factor,subdiv_factor); % create sub image
        %so that # of subpixels that can potentially be labeled is:
        n_to_label=size(w_a,1)*subdiv_factor^3;
        
        if n_to_label<n_mol_per_cluster
            disp(['Fluorophore: ',sm_par.fluorophore_name]);
            disp(['Number of possibly labeled pixels: ',num2str(n_to_label)]);
            errordlg('Number of molecules to label is too high; reduce # of molecules','SMIS Error');
            uiwait
            ok=0;
            return
        end
        
        %pick randomly locations where labeling molecules will end up
        r=rand(1,n_to_label); %uniformly distributed random number
        [~, sorted_i]=sort(r,'ascend'); %sorted_i are the indices of sorted_r
        sorted_i=sorted_i'; %make it a line vector
        
        % below: corrected bug 22/04/2021
        sorted_i_a = fix((sorted_i-1)/subdiv_factor^3); %the index of the big pixel in  a
        sorted_i_r = fix(sorted_i-sorted_i_a*subdiv_factor^3) ; %the index of the small pixel in sub_a
        sorted_i_a=sorted_i_a(1:n_mol_per_cluster)+1;
        sorted_i_r=sorted_i_r(1:n_mol_per_cluster);
        if max(sorted_i_a) > length(w_a)
            disp(['Fluorophore: ',sm_par.fluorophore_name]);
            errordlg('Retry or reduce nearest neighbour distance (try 0.25 nm !) !','SMIS Error');
            uiwait
            ok=0;
            return
        end
        
        i_selected=w_a(sorted_i_a);
        for j=1:n_mol_per_cluster
            a_sm(i_selected(j))=a_sm(i_selected(j))+100; %Loop to figure out if multiple SM in 1 pixel
        end
        [x0,y0,z0]=ind2sub(size(a),i_selected); % x is the row #; y is the column #
        [sub_x,sub_y,sub_z]=ind2sub(size(sub_a),sorted_i_r);
        
        % Get the xy position of all molecules and add some random components
        new_x=x0-0.5+1/subdiv_factor*((sub_x-0.5)+(rand(n_mol_per_cluster,1)-0.5)); % molecule in pixel x might take positions from x-0.5 to x+0.5
        new_y=y0-0.5+1/subdiv_factor*((sub_y-0.5)+(rand(n_mol_per_cluster,1)-0.5)); % molecule in pixel x might take positions from y-0.5 to y+0.5
        new_z=z0-0.5+1/subdiv_factor*((sub_z-0.5)+(rand(n_mol_per_cluster,1)-0.5)); % molecule in pixel x might take positions from y-0.5 to y+0.5
        x=[x; new_x];
        y=[y; new_y];
        z=[z; new_z];
    end
    
    % Place the rest of the molecules out of the pattern
    %indices of pixels that might be labeled
    if n_mol_op > 0
        w_a=find(a==0);
        [x_op,y_op,z_op,a_sm] = place_single_molecules_on_pattern(a_sm,w_a, im_par,n_mol_op);
        x=[x; x_op];
        y=[y; y_op];
        z=[z; z_op];
    end
end
end