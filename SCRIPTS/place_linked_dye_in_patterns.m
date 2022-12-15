function [sm_ref, sm_target, sm_par_target, a_target_sm, n_mol_ref_out, n_mol_target_out,ok] = place_linked_dye_in_patterns(sm_ref, ...
    sm_target, sm_par_target, a_ref, a_target, im_par, n_mol_ref, n_mol_target, pair_id)
% NAME:
%	PLACE_LINKED_DYES_IN_PATTERNS
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
%   sm_par_target: sm_par of target
%	a_ref: the pattern image onto place the linked target SM
%	a_target: the pattern image onto place the non paired target SM
%	im_par: some general parameters
%   n_mol_ref: # of ref SM 
%   n_mol_target: # of target SM 
%   pair_id: the index of fluorophore_pairs corresponding to the reference and target SMs (index in sm_par)
%
% OUTPUTS:
%   sm_target: the updated target SM
%   sm_par_target: updated sm_par of target
%   a_target_sm = the image, where the target SMs have been positionned
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
%	D.Bourgeois, November 2019: version > simulate_palm_vsn15.3
%	D.Bourgeois, March 2021: Added error handling

ok=1;
s=min([n_mol_ref, n_mol_target]);

r1=im_par.d1d2_dist(pair_id)*im_par.binning/im_par.raster; %in pixels for the unbinned image
r2=im_par.d1d2_dist_sig(pair_id)*im_par.binning/im_par.raster; %in pixels for the unbinned image
constrained=im_par.d2_constrained(pair_id);
x2=zeros(n_mol_target,1);
y2=zeros(n_mol_target,1);
c_sp_target=zeros(n_mol_target,1);
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
            x2(k) = x(k) + d.*cos(theta);
            y2(k) = y(k) + d.*sin(theta);
            %Check if new coordinates are inside a_2 image
            if constrained
                if ~(x2(k) < 1 || y2(k) < 1 || x2(k) > size(a_ref,1) || y2(k) > size(a_ref,2))     
                    if a_ref(round(x2(k)),round(y2(k)))==c_sp(k)
                        placing_ok=1;
                        c_sp_target(k)=c_sp(k); % update subpattern id of target
                    else
                        x2(k)=1e+20;
                        y2(k)=1e+20;
                    end
                end
            else
                placing_ok=1;
            end
            j=j+1;
        end
        
        if ~placing_ok % Failure
            x2(k)=1e+20;
            y2(k)=1e+20;
        end
    end
       
    if n_mol_target > n_mol_ref % Eventually place additionnal target molecules on target patterns
        [x2(n_mol_ref+1:n_mol_target),y2(n_mol_ref+1:n_mol_target),~,...
            c_sp_target(n_mol_ref+1:n_mol_target),~,sm_par_target,ok] = ...
            place_additionnal_single_molecules_simple(a_target,im_par,sm_par_target, n_mol_target-n_mol_ref);
        if ~ok; return; end
    end
    
    %Initialize a_target_sm
    a_target_sm=a_target;
    a_target_sm(:)=0;

    %remove target molecules out of the field of view
    w_out=find(x2 < 0.5 | x2 >= (size(a_target,1)+0.5) | y2 < 0.5 | y2 >= (size(a_target,2)+0.5));
    n_mol_target_out=n_mol_target-length(w_out); % update # of target molecules
    if n_mol_target_out==0
        errordlg('No target molecules could be placed !','SMIS Error');
        uiwait
        ok=0;
        return
    end

    sm_target(w_out)=[];
    x2(w_out)=[]; y2(w_out)=[];
    c_sp_target(w_out)=[]; 
    newVals = num2cell(x2); [sm_target.x] = newVals{:}; 
    newVals = num2cell(y2); [sm_target.y] = newVals{:};
    newVals = num2cell(c_sp_target); [sm_target.c_sp] = newVals{:};
            
    %update sm_ref
    w_out_ref=w_out(w_out<=n_mol_ref);
    sm_ref(w_out_ref)=[];
    n_mol_ref_out=length(sm_ref); % update # of ref molecules
    
    removed = n_mol_target-n_mol_target_out; % # of removed molecules
    
    if removed>0
        placing_ok = questdlg(['WARNING: ',num2str(removed), ' target molecules could not be placed ! Do you want to continue ? '], 'SMIS Question','Yes','No','Yes');
        if strcmp(placing_ok,'No')==1
            [sm_target.x]=deal(nan);
            [sm_target.y]=deal(nan);
            [sm_target.z]=deal(nan);
            [sm_target.c_sp]=deal(nan);
            [sm_target.matched]=deal(nan);
            a_target_sm=nan;
            ok=0;
            return
        end
    end
    
    for i=1:n_mol_target_out
        %Loop to figure out if multiple SM in 1 pixel
        a_target_sm(round(sm_target(i).x),round(sm_target(i).y))=a_target_sm(round(sm_target(i).x),round(sm_target(i).y))+100;
    end
    
%go to 3D case    
elseif im_par.simul_3D==1 % 3D mode
    z2=zeros(n_mol_target,1);
    z=[sm_ref.z];
    max_trials=1000;
    for k=1:s
        placing_ok=0;
        j=1;
        while ~placing_ok && j<max_trials
            d = r1 + r2*randn(1,1); % define random distances
            theta = 2*pi*rand(1,1); % define random angle
            phi = -pi/2+pi*rand(1,1); % define random angle
            x2(k) = x(k) + d*cos(theta)*cos(phi);
            y2(k) = y(k) + d*sin(theta)*cos(phi);
            z2(k) = z(k) + d*sin(phi);
            %Check if new coordinates are inside a_2 image
            if constrained
                if ~(x2(k) < 1 || y2(k) < 1 || z2(k) < 1 || x2(k) > size(a_ref,1) || y2(k) > size(a_ref,2) || z2(k) > size(a_ref,3))
                    if a_ref(round(x2(k)),round(y2(k)),round(z2(k)))==c_sp(k)
                        placing_ok=1;
                        c_sp_target(k)=c_sp(k); % update subpattern id of target
                    else
                        x2(k)=1e+20;
                        y2(k)=1e+20;
                        z2(k)=1e+20;
                    end
                end
                j=j+1;
            else
                placing_ok=1;
            end
            
        end
        
        if ~placing_ok % Failure
            x2(k)=1e+20;
            y2(k)=1e+20;
            z2(k)=1e+20;
        end

    end
    
    if n_mol_target > n_mol_ref % Eventually Place additionnal molecules
        [x2(n_mol_ref+1:n_mol_target),y2(n_mol_ref+1:n_mol_target),~,...
            c_sp_target(n_mol_ref+1:n_mol_target),~,sm_par_target,ok] = ...
            place_additionnal_single_molecules_simple(a_target,im_par,sm_par_target, n_mol_target-n_mol_ref);
        if ~ok; return; end
    end
    
    %Initialize a_target_sm
    a_target_sm=a_target;
    a_target_sm(:)=0;
    
    %remove molecules out of the field of view
    w_out=find(x2 < 0.5 | x2 >= (size(a_target,1)+0.5) | y2 < 0.5 | y2 >= (size(a_target,2)+0.5) | z2 < 0.5 | z2 >= (size(a_target,3)+0.5));
    n_mol_target_out=n_mol_target-length(w_out); % update # of target molecules

    if n_mol_target_out==0
        errordlg('No target molecules could be placed !','SMIS Error');
        uiwait
        ok=0;
        return
    end

    sm_target(w_out)=[];
    x2(w_out)=[]; y2(w_out)=[]; z2(w_out)=[];
    c_sp_target(w_out)=[]; 
    newVals = num2cell(x2); [sm_target.x] = newVals{:}; 
    newVals = num2cell(y2); [sm_target.y] = newVals{:};
    newVals = num2cell(z2); [sm_target.z] = newVals{:};
    newVals = num2cell(c_sp_target); [sm_target.c_sp] = newVals{:};
            
    %update sm_ref
    w_out_ref=w_out(w_out<=n_mol_ref);
    sm_ref(w_out_ref)=[];
    n_mol_ref_out=length(sm_ref); % update # of ref molecules
   
    removed = n_mol_target-n_mol_target_out; % # of removed molecules
    
    if removed>0
        placing_ok = questdlg(['WARNING: ',num2str(removed), ' molecules could not be placed ! Do you want to continue ? '], 'SMIS Question','Yes','No','Yes');
        if strcmp(placing_ok,'No')==1
            [sm_target.x]=deal(nan);
            [sm_target.y]=deal(nan);
            [sm_target.z]=deal(nan);
            [sm_target.c_sp]=deal(nan);
            [sm_target.matched]=deal(nan);
            a_target_sm=nan;
            ok=0;
            return
        end
    end
    
    for i=1:(n_mol_target-removed)
        %Loop to figure out if multiple SM in 1 pixel
        a_target_sm(round(sm_target(i).x),round(sm_target(i).y),round(sm_target(i).z))=a_target_sm(round(sm_target(i).x),round(sm_target(i).y),round(sm_target(i).z))+200;       
    end
end
end