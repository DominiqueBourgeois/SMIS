function [x,y,z,a_sm,ok] = add_single_molecules_on_pattern(a_sm, w_a,sm_par, im_par,n_mol)

% PURPOSE:
%	Add molecules on a_sm image on a specified pattern
%
% CATEGORY:
%	Single molecule.
%
% INPUTS:
%	a_sm: the image onto place the SMs
%   w_a: the indices where to place the molecules
%   sm_par: single molecule parameters
%	im_par: imaging parameters
%   n_mol: the # of molecules to add
%
% OUTPUTS:
%	x = the x coordinates [in raster units] on high-resolution image
%   y = the y coordinates [in raster units] on high-resolution image
%   z = the z coordinates [in raster units] on high-resolution image
%   a_sm = the updated image, where the SMs have been positionned
%   ok = 1 if execution was okay, 0 if an error occurred
%
% SIDE EFFECTS:
%	None.
%
% RESTRICTIONS:
%	None.
%
% MODIFICATION HISTORY:
%	D.Bourgeois, November 2019:
%	D.Bourgeois, March 2021: Added error handling

ok=1;
max_n_to_label=1e+08;

if numel(size(a_sm))==2 % 2D image
    
    %each pixel should be subdivided to account for nearest neighbour distance
    subdiv_factor=max([1,fix(im_par.raster/sm_par.near_dist/im_par.binning)]);
    sub_a=zeros(subdiv_factor,subdiv_factor); % create sub image
    %so that # of subpixels that can potentially be labeled is:
    n_to_label=size(w_a,1)*subdiv_factor^2;
    
    if n_to_label<n_mol
        disp(['Fluorophore: ',sm_par.fluorophore_name]);
        disp(['Number of possibly labeled pixels for subpattern: ',num2str(n_to_label)]);
        errordlg('Number of molecules to label is too high; reduce # of molecules','SMIS Error');
        uiwait
        ok=0;
        return
    end
    
    if n_to_label>max_n_to_label
        disp(['Fluorophore: ',sm_par.fluorophore_name]);
        disp(['Number of possibly labeled pixels for subpattern is too high for memory: ',num2str(n_to_label)]);
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
    %sorted_i_a = fix(sorted_i/subdiv_factor^2); %the index of the big pixel in  a
    %sorted_i_r = fix(sorted_i-sorted_i_a*subdiv_factor^2) ; %the index of the small pixel in sub_a
    %sorted_i_a=sorted_i_a(1:n_mol)+1;
    %sorted_i_r=sorted_i_r(1:n_mol)+1;
    if max(sorted_i_a) > length(w_a)
        disp('***** Problem with dense labeling for subpattern ********');
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
    x=x0-0.5+1/subdiv_factor*((sub_x-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel x might take positions from x-0.5 to x+0.5
    y=y0-0.5+1/subdiv_factor*((sub_y-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel x might take positions from y-0.5 to y+0.5
    z=0*x+1; % Set z to 1 for 2D
    
elseif numel(size(a))==3 % 3D image
    
    %each pixel should be subdivided to account for nearest neighbour distance
    subdiv_factor=max([1,fix(im_par.raster/sm_par.near_dist/im_par.binning)]);
    sub_a=zeros(subdiv_factor,subdiv_factor,subdiv_factor); % create 3D sub image
    %so that # of subpixels that can potentially be labeled is:
    n_to_label=size(w_a,1)*subdiv_factor^3;
    
    if n_to_label<n_mol
        disp(['Fluorophore: ',sm_par.fluorophore_name]);
        disp(['Number of possibly labeled pixels for subpattern: ',num2str(n_to_label)]);
        error('Number of molecules to label is too high; reduce # of molecules');
    end
    
    if n_to_label>max_n_to_label
        disp(['Fluorophore: ',sm_par.fluorophore_name]);
        disp(['Number of possibly labeled pixels for subpattern is too high for memory: ',num2str(n_to_label)]);
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
    %sorted_i_a = fix(sorted_i/subdiv_factor^3); %the index of the big pixel in  a
    %sorted_i_r = fix(sorted_i-sorted_i_a*subdiv_factor^3) ; %the index of the small pixel in sub_a
    %sorted_i_a=sorted_i_a(1:n_mol)+1;
    %sorted_i_r=sorted_i_r(1:n_mol)+1;    
    if max(sorted_i_a) > length(w_a)
        disp('***** Problem with dense labeling for subpattern ********');
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
    x=x0-0.5+1/subdiv_factor*((sub_x-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel x might take positions from x-0.5 to x+0.5
    y=y0-0.5+1/subdiv_factor*((sub_y-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel y might take positions from y-0.5 to y+0.5
    z=z0-0.5+1/subdiv_factor*((sub_z-0.5)+(rand(n_mol,1)-0.5)); % molecule in pixel z might take positions from z-0.5 to z+0.5
end

