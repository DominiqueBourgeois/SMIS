function [sm_ref, sm_target, a_sm, n_mol_ref_out, n_mol_target_out,ok] = place_linked_dye_in_two_patterns(sm_ref, ...
    sm_target, a_ref, im_par, n_mol_ref, n_mol_target, pair_id)
% NAME:
%	PLACE_LINKED_DYES
%
% PURPOSE:
%	For a covalently bound SM pair, get the xyz coordinates of the Target
%	SM from knowledge of the positions of ref SM
%
% CATEGORY:
%	Single molecule.
%
% INPUTS:
%   sm_ref: the reference SM
%   sm_target: the target SM
%	a_ref: the image onto place the target SM
%	im_par: some general parameters
%   n_mol_ref: # of ref SM 
%   n_mol_target: # of target SM 
%   pair_id: the index of fluorophore_pairs corresponding to the reference and target SMs (index in sm_par)
%
% OUTPUTS:
%   sm_target: the updated target SM
%   a_sm = the image, where the target SMs have been positionned
%   n_mol_target: updated # of target SM 
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
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, March 2021: Added error handling

ok=1;
s=min([n_mol_ref, n_mol_target]);

r1=im_par.d1d2_dist(pair_id)*im_par.binning/im_par.raster; %in pixels for the unbinned image
r2=im_par.d1d2_dist_sig(pair_id)*im_par.binning/im_par.raster; %in pixels for the unbinned image
constrained=im_par.d2_constrained(pair_id);
x_2=zeros(n_mol_target,1);
y_2=zeros(n_mol_target,1);
x=[sm_ref.x];
y=[sm_ref.y];
c_sp=[sm_ref.c_sp];

if im_par.simul_3D==0 % 2D mode
    max_trials=1000;
    for k=1:s
        placing_ok=0;
        j=1;
        while ~placing_ok && j<max_trials
            d = r1 + r2*randn(1,1); % define random distances
            theta = 2*pi*rand(1,1); % define random angle
            x_2(k) = x(k) + d.*cos(theta);
            y_2(k) = y(k) + d.*sin(theta);
            %Check if new coordinates are inside a_2 image
            if constrained
                if ~(x_2(k) < 1 || y_2(k) < 1 || x_2(k) > size(a_ref,1) || y_2(k) > size(a_ref,2))     
                    if a_ref(round(x_2(k)),round(y_2(k)))==c_sp(k)
                        placing_ok=1;
                    else
                        x_2(k)=1e+20;
                        y_2(k)=1e+20;
                    end
                end
            else
                placing_ok=1;
            end
            j=j+1;
        end
        
        if ~placing_ok % Failure
            x_2(k)=1e+20;
            y_2(k)=1e+20;
        end
    end
       
    if n_mol_target > n_mol_ref % Eventually place additionnal target molecules
        [x_2(n_mol_ref+1:n_mol_target),y_2(n_mol_ref+1:n_mol_target),~, ~] = ...
            place_single_molecules_simple(a_ref,im_par,n_mol_target-n_mol_ref);
    end
    
    a_sm=a_ref;
    a_sm(:)=0;
    
    %remove target molecules out of the field of view
    w_ok=find(x_2 >= 0.5 & x_2 < (size(a_ref,1)+0.5) & y_2 >= 0.5 & y_2 < (size(a_ref,2)+0.5)); 
    n_mol_target_out=length(w_ok); % update # of target molecules
    if n_mol_target_out==0
        errordlg('No target molecules could be placed !','SMIS Error');
        uiwait
        ok=0;
        return
    end

    sm_target=sm_target(w_ok);
    newVals = num2cell(x_2(w_ok)); [sm_target.x] = newVals{:}; 
    newVals = num2cell(y_2(w_ok)); [sm_target.y] = newVals{:};
            
    %update sm_ref
    w_ok_ref=w_ok(w_ok<=n_mol_ref);
    n_mol_ref_out=length(w_ok_ref); % update # of target molecules
    sm_ref=sm_ref(w_ok_ref);
    
    removed = n_mol_target-n_mol_target_out; % # of removed molecules
    
    if removed>0
        placing_ok = questdlg(['WARNING: ',num2str(removed), ' molecules could not be placed ! Do you want to continue ? '], 'SMIS Question','Yes','No','Yes');
        if strcmp(placing_ok,'No')==1
            [sm_target.x]=deal(nan);
            [sm_target.y]=deal(nan);
            [sm_target.z]=deal(nan);
            [sm_target.matched]=deal(nan);
            a_sm=nan;
            ok=0;
            return
        end
    end
   
    for i=1:n_mol_target_out
        %Loop to figure out if multiple SM in 1 pixel
        a_sm(round(sm_target(i).x),round(sm_target(i).y))=a_sm(round(sm_target(i).x),round(sm_target(i).y))+200;
    end
elseif im_par.simul_3D==1 % 3D mode
    z_2=zeros(n_mol_target,1);
    z=[sm_ref.z];
    max_trials=1000;
    for k=1:s
        placing_ok=0;
        j=1;
        while ~placing_ok && j<max_trials
            d = r1 + r2*randn(1,1); % define random distances
            theta = 2*pi*rand(1,1); % define random angle
            phi = -pi/2+pi*rand(1,1); % define random angle
            x_2(k) = x(k) + d*cos(theta)*cos(phi);
            y_2(k) = y(k) + d*sin(theta)*cos(phi);
            z_2(k) = z(k) + d*sin(phi);
            %Check if new coordinates are inside a_2 image
            if constrained
                if ~(x_2(k) < 1 || y_2(k) < 1 || z_2(k) < 1 || x_2(k) > size(a_ref,1) || y_2(k) > size(a_ref,2) || z_2(k) > size(a_ref,3))
                    if a_ref(round(x_2(k)),round(y_2(k)),round(z_2(k)))>0
                        placing_ok=1;
                    else
                        x_2(k)=1e+20;
                        y_2(k)=1e+20;
                        z_2(k)=1e+20;
                    end
                end
                j=j+1;
            else
                placing_ok=1;
            end
            
        end
        
        if ~placing_ok % Failure
            x_2(k)=1e+20;
            y_2(k)=1e+20;
            z_2(k)=1e+20;
        end

    end
    
    if n_mol_target > n_mol_ref % Eventually Place additionnal molecules
        [x_2(n_mol_ref+1:n_mol_target),y_2(n_mol_ref+1:n_mol_target),z_2(n_mol_ref+1:n_mol_target), ~] = ...
            place_single_molecules_simple(a_ref,im_par,n_mol_target-n_mol_ref);
    end
    
    a_sm=a_ref;
    a_sm(:)=0;
    
    %remove molecules out of the field of view
    w_ok=find(x_2 >= 0.5 & x_2 < (size(a_ref,1)+0.5) & y_2 >= 0.5 & y_2 < (size(a_ref,2)+0.5) & z_2 >= 0.5 & z_2 < (size(a_ref,3)+0.5));
    n_mol_target_out=length(w_ok); % update # of target molecules
    if n_mol_target_out==0
        errordlg('No target molecules could be placed !','SMIS Error');
        uiwait
        ok=0;
        return
    end

    sm_target=sm_target(w_ok);
    newVals = num2cell(x_2(w_ok)); [sm_target.x] = newVals{:};
    newVals = num2cell(y_2(w_ok)); [sm_target.y] = newVals{:};
    newVals = num2cell(z_2(w_ok)); [sm_target.z] = newVals{:};
        
    %update sm_ref
    w_ok_ref=w_ok(w_ok<=n_mol_ref);
    n_mol_ref_out=length(w_ok_ref); % update # of target molecules   
    sm_ref=sm_ref(w_ok_ref);

    removed = n_mol_target-n_mol_target_out; % # of removed molecules
    
    if removed>0
        placing_ok = questdlg(['WARNING: ',num2str(removed), ' molecules could not be placed ! Do you want to continue ? '], 'SMIS Question','Yes','No','Yes');
        if strcmp(placing_ok,'No')==1
            [sm_target.x]=deal(nan);
            [sm_target.y]=deal(nan);
            [sm_target.z]=deal(nan);
            [sm_target.matched]=deal(nan);
            a_sm=nan;
            ok=0;
            return
        end
    end
    
    for i=1:(n_mol_target-removed)
        %Loop to figure out if multiple SM in 1 pixel
        a_sm(round(sm_target(i).x),round(sm_target(i).y),round(sm_target(i).z))=a_sm(round(sm_target(i).x),round(sm_target(i).y),round(sm_target(i).z))+200;       
    end
end
end