function [x,y,z,sp,a_sm,sm_par,ok] = place_additionnal_single_molecules_simple(a,im_par,sm_par,add_n_mol)
% NAME:
%	PLACE_ADDITIONAL_SINGLE_MOLECULES
%
% PURPOSE:
%	Get xyz coordinates of additional SMs following a pattern on an image
%
% CATEGORY:
%	Single molecule.
%
% CALLING SEQUENCE:
%	[x,y,z,sp,a_sm,sm_par] = place_single_molecules(a,im_par,sm_par,n_mol)
%
% INPUTS:
%	a: the image onto place the SMs
%	im_par: imaging parameters
%   sm_par: single molecule parameters
%   n_mol: the additional # of molecules to place
%
% OUTPUTS:
%	x = the x coordinates [in raster units] on high-resolution image
%   y = the y coordinates [in raster units] on high-resolution image
%   z = the z coordinates [in raster units] on high-resolution image
%   sp = the subpattern # on which each SM is placed
%   a_sm = the image, where the SMs have been positionned
%   sm_par = updated sm_par
%   ok = 1 if execution was okay, 0 if an error occurred
%
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
%	D.Bourgeois, November 2019: place molecules according to subpatterns
%	D.Bourgeois, March 2021: Added error handling

ok=1;
max_n_to_label=1e+08;

n_sp_id=sm_par.n_sp_id;
n_sp=numel(n_sp_id); % # of subpatterns

sm_par.w_patterns(1:n_sp)=struct('w', []); % will contain the indices of subpatterns in image. 

%define x, y, and z
x=zeros(add_n_mol,1);
y=zeros(add_n_mol,1);
z=zeros(add_n_mol,1);
sp=zeros(add_n_mol,1); 

a_sm=0*a; % Set a_sm

mol_id=0; % initiate mol counter

% Fraction of molecules in different subpatterns, the first corresponds to
% molecules out of any pattern
fr=sm_par.n_sp_fr;

for k=1:n_sp
    n_mol=fix(add_n_mol*fr(k)); % # of molecules to be placed on that subpattern

    %indices of pixels that might be labeled
    w_a=find(a==n_sp_id(k));
    sm_par.w_patterns(k).w=w_a; % update sm_par
    
    if n_mol>0   
        if numel(size(a))==2 % 2D image
            
            %each pixel should be subdivided to account for nearest neighbour distance
            subdiv_factor=max([1,fix(im_par.raster/sm_par.near_dist/im_par.binning)]);
            sub_a=zeros(subdiv_factor,subdiv_factor); % create sub image
            %so that # of subpixels that can potentially be labeled is:
            n_to_label=size(w_a,1)*subdiv_factor^2;
            
            if n_to_label<n_mol
                disp(['Fluorophore: ',sm_par.fluorophore_name]);
                disp(['Number of possibly labeled pixels for subpattern ', num2str(n_sp_id(k)) ,' : ',num2str(n_to_label)]);
                errordlg('Number of molecules to label is too high; reduce # of molecules','SMIS Error');
                uiwait
                ok=0;
                return
            end
            
            if n_to_label>max_n_to_label
                disp(['Fluorophore: ',sm_par.fluorophore_name]);
                disp(['Number of possibly labeled pixels for subpattern ', num2str(n_sp_id(k)), ' is too high for memory: ',num2str(n_to_label)]);
                errordlg('Increase ''near_dist'', or decrease ''binning'' or ''raster'' parameters','SMIS Error');
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
            sorted_i_a=sorted_i_a(1:n_mol)+1;
            sorted_i_r=sorted_i_r(1:n_mol);
%             sorted_i_a = fix(sorted_i/subdiv_factor^2); %the index of the big pixel in  a
%             sorted_i_r = fix(sorted_i-sorted_i_a*subdiv_factor^2) ; %the index of the small pixel in sub_a
%             sorted_i_a=sorted_i_a(1:n_mol)+1;
%             sorted_i_r=sorted_i_r(1:n_mol)+1;
            if max(sorted_i_a) > length(w_a)
                disp(['***** Problem with dense labeling for subpattern ', num2str(n_sp_id(k)), '********']);
                disp(['Fluorophore: ',sm_par.fluorophore_name]);
                errordlg('Retry or reduce nearest neighbour distance (should be smaller than raster size/2)!','SMIS Error');
                uiwait
                ok=0;
                return
            end
            
            i_selected=w_a(sorted_i_a);
            for i=1:n_mol
                a_sm(i_selected(i))=a_sm(i_selected(i))+100; %Loop to figure out if multiple SM in 1 pixel
            end
            [x0,y0]=wheresub(i_selected,a); % x is the row #; y is the column #
            [sub_x,sub_y]=wheresub(sorted_i_r,sub_a);
            % Get the xy position of all molecules and add some random components
            x(mol_id+1:mol_id+n_mol)=x0-0.5+1/subdiv_factor*((sub_x-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel x might take positions from x-0.5 to x+0.5
            y(mol_id+1:mol_id+n_mol)=y0-0.5+1/subdiv_factor*((sub_y-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel x might take positions from y-0.5 to y+0.5
            z(mol_id+1:mol_id+n_mol)=0*x(mol_id+1:mol_id+n_mol)+1; % Set z to 1 for 2D
            sp(mol_id+1:mol_id+n_mol)=n_sp_id(k); % Index of the subpattern
        elseif numel(size(a))==3 % 3D image
            
            %each pixel should be subdivided to account for nearest neighbour distance
            subdiv_factor=max([1,fix(im_par.raster/sm_par.near_dist/im_par.binning)]);
            sub_a=zeros(subdiv_factor,subdiv_factor,subdiv_factor); % create 3D sub image
            %so that # of subpixels that can potentially be labeled is:
            n_to_label=size(w_a,1)*subdiv_factor^3;
            
            if n_to_label<n_mol
                disp(['Fluorophore: ',sm_par.fluorophore_name]);
                disp(['Number of possibly labeled pixels for subpattern ', num2str(n_sp_id(k)) ,' : ',num2str(n_to_label)]);
                errordlg('Number of molecules to label is too high; reduce # of molecules','SMIS Error');
                uiwait
                ok=0;
                return
            end
            
            if n_to_label>max_n_to_label
                disp(['Fluorophore: ',sm_par.fluorophore_name]);
                disp(['Number of possibly labeled pixels for subpattern ', num2str(n_sp_id(k)), ' is too high for memory: ',num2str(n_to_label)]);
                errordlg('Increase ''near_dist'', or decrease ''binning'' or ''raster'' parameters','SMIS Error');
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
            sorted_i_r = fix(sorted_i-sorted_i_a*subdiv_factor^3); %the index of the small pixel in sub_a
            sorted_i_a=sorted_i_a(1:n_mol)+1;
            sorted_i_r=sorted_i_r(1:n_mol);
%             sorted_i_a = fix(sorted_i/subdiv_factor^3); %the index of the big pixel in  a
%             sorted_i_r = fix(sorted_i-sorted_i_a*subdiv_factor^3) ; %the index of the small pixel in sub_a
%             sorted_i_a=sorted_i_a(1:n_mol)+1;
%             sorted_i_r=sorted_i_r(1:n_mol)+1;
            if max(sorted_i_a) > length(w_a)
                disp(['***** Problem with dense labeling for subpattern ', num2str(n_sp_id(k)), '********']);
                disp(['Fluorophore: ',sm_par.fluorophore_name]);
                errordlg('Retry or reduce nearest neighbour distance (should be smaller than raster size/2)!','SMIS Error');
                uiwait
                ok=0;
                return
            end
            
            i_selected=w_a(sorted_i_a);
            for i=1:n_mol
                a_sm(i_selected(i))=a_sm(i_selected(i))+100; %Loop to figure out if multiple SM in 1 pixel
            end
            [x0,y0,z0]=ind2sub(size(a),i_selected); % x is the row #; y is the column #
            [sub_x,sub_y,sub_z]=ind2sub(size(sub_a),sorted_i_r);
            % Get the xyz position of all molecules and add some random components
            x(mol_id+1:mol_id+n_mol)=x0-0.5+1/subdiv_factor*((sub_x-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel x might take positions from x-0.5 to x+0.5
            y(mol_id+1:mol_id+n_mol)=y0-0.5+1/subdiv_factor*((sub_y-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel y might take positions from y-0.5 to y+0.5
            z(mol_id+1:mol_id+n_mol)=z0-0.5+1/subdiv_factor*((sub_z-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel z might take positions from z-0.5 to z+0.5
            sp(mol_id+1:mol_id+n_mol)=n_sp_id(k); % Index of the subpattern
        end
        mol_id=mol_id+n_mol; % increment molecule counter
    end
end