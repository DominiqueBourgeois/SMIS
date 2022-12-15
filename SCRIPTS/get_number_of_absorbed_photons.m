function N=get_number_of_absorbed_photons(sm, lasers, sm_par, im_par, sampling_rate)

% PURPOSE:
%   Get the number of absorbed photon per sampling time from each photoactive state according
%   to single molecule position, extinction coefficient, laser power
%   density 
%
% INPUTS:
%   sm: the single molecule (with coordinates on high-resolution image) in raster units
%   lasers: the lasers
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   N: number of photons absorbed per sampling time 
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, June 2020. Changed from # of photons per second to # of
%	photons per sampling time
%	D.Bourgeois, November 2020. Added TIRF or HILO modes
%   D.Bourgeois, November 2021, Corrected a bug for N<1 photons: the use of "round" was incorrect !!

avogadro = 6.02e+23; % [mol-1]

n_lasers=size(lasers,2);
update_lasers=0;

%Check if we need to recalculate
if im_par.add_drift==0 && im_par.add_diffusion==0 && im_par.current_frame>1
    for i=1:n_lasers
        if lasers(i).sequence(im_par.current_frame)~=lasers(i).sequence(im_par.current_frame-1)
            update_lasers=1; % Then update Ns
        end
    end
    if update_lasers==0 % Then do not update Ns
        if im_par.during_frametime==1
            N=sm.Ns(end,:);
        else
            N=sm.Ns(1,:);
        end
        return
    end
end

% Get sm coordinates
if im_par.simul_3D==1
    [x,y,z]=get_coordinates_on_detector(sm, im_par); % x,y,z in raster
    % do nothing if the molecule is out of the FOV including in Z (for
    % diffusion measurements, care should be taken to choose a sufficiently thick sample 
    if round(x)<=0 || round(x)>im_par.n || round(y)<=0 || round(y)>im_par.m || round(z)<=0 || round(z)>im_par.nz % only calculate if molecule within FOV
        N=zeros(1,sm_par.n_states); % No photon absorbed
        return
    end
else
    [x,y,~]=get_coordinates_on_detector(sm, im_par); % x,y,z in raster
    % do nothing if the molecule is out of the FOV
    if round(x)<=0 || round(x)>im_par.n || round(y)<=0 || round(y)>im_par.m % only calculate if molecule within FOV
        N=zeros(1,sm_par.n_states); % No photon absorbed
        return
    end
end

%Define the array of absorbed photons
N=zeros(n_lasers,sm_par.n_states); % 

%Go over all lasers
for i=1:n_lasers
    if (im_par.during_frametime==0 && lasers(i).on_during_addtime==1) || (im_par.during_frametime==1 && lasers(i).on_during_frametime==1)
        if lasers(i).sequence(im_par.current_frame) > 0 && ~isempty(lasers(i).beam_profile)
            %Local laser energy density [photons/cm2/s]           
            local_density = lasers(i).sequence(im_par.current_frame)/100*...
                lasers(i).beam_profile(round(x),round(y))*(1e+7/im_par.raster)^2; % 1cm=1e+7nm

            %Treat TIRF or HILO cases in 3D mode
            if im_par.simul_3D==1 && lasers(i).tirf==1
                local_density=local_density*get_tirf_hilo_correction_factor(z, lasers(i), im_par);
            end
                       
            %Go over fluorescent states
            for j=1:sm_par.n_fluorescent_states
                fluo_state_id=sm_par.fluorescent_states(j);
                
                % Only scan for possible states if molecule already photoconverted
                if sm.state<sm_par.converted_state || (sm.state>=sm_par.converted_state && fluo_state_id >= sm_par.converted_state)
                    
                    %The cross section is given by
                    eps=sm_par.spectral_data.exc_spectra(j).eps_lasers(i);
                    cross_section = 1000*log(10)*eps/avogadro; % [cm2]

                    if im_par.apply_poisson_stat_for_lasers==1
                        N(i,fluo_state_id)=poissrnd(local_density*cross_section); % 
                    else
                        %N(i,fluo_state_id)=round(local_density*cross_section); % Less accurate but more rapid. Wrong for N<1 !!
                        N(i,fluo_state_id)=local_density*cross_section; % Less accurate but more rapid. 
                    end                    
                end
            end
            
            %go over photoactive dark states
            for j=1:sm_par.n_photoactive_dark_states
                dark_state_id=sm_par.photoactive_dark_states(j);
                
                % Only scan for possible states if molecule already photoconverted
                if sm.state<sm_par.converted_state || (sm.state>=sm_par.converted_state && dark_state_id >= sm_par.converted_state)
                    %The cross section is given by
                    eps=sm_par.spectral_data.dark_spectra(j).eps_lasers(i);
                    cross_section = 1000*log(10)*eps/avogadro; % [cm2]                  
                    
                    if im_par.apply_poisson_stat_for_lasers==1
                        N(i,dark_state_id)=poissrnd(local_density*cross_section); % 
                    else
                        %N(i,dark_state_id)=round(local_density*cross_section); % Less accurate but more rapid. Wrong for N<1 !!
                        N(i,dark_state_id)=local_density*cross_section; % Less accurate but more rapid.
                    end
                end
            end
            
            %do the correction for anisotropy
            if sm_par.anisotropy==1 % take anisotropy into account
                theta=sm.theta; phi=sm.phi;
                if lasers(i).polarization==1 % laser is linearly polarized
                    %the factor 3 comes from the fact that
                    %<cos((phi-act_laser_par.phi_laser)*pi/180)^2>=1/2
                    %and that
                    %<cos(theta*pi/180)^2>=2/3
                    %so 1/2*2/3=1/3 so for randomly oriented molecules we have to get
                    %no modification of cross_section, hence the factor of 3
                    N(i,:)=N(i,:)*2*cos((phi-lasers(i).phi)*pi/180)^2*3/2*cos(theta*pi/180)^2;
                else % laser is circularly polarized
                    %here the factor 3/2 comes from the fact that
                    %<cos(theta*pi/180)^2>=2/3
                    N(i,:)=N(i,:)*3/2*cos(theta*pi/180)^2;
                end
            end
        end
    end
end

%Finally add up all absorbed photons by the various lasers and get N per
%sampling time
N=(1/sampling_rate)*sum(N,1);

