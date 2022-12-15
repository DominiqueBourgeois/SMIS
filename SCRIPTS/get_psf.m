function sm_par=get_psf(n_fluorophores, sm_par, im_par, psf_mode, channel, simul_3D)
% NAME:
%	get_psf
%
% PURPOSE:
%       Simulates a PSF for each fluorophore
%
% INPUTS:
%   n_fluorophores: the # of fluorophores
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%   psf_mode: the type of algorithm to get the psf (only Gaussian for now)
%   channel: 1 or 2
%   simul_3D: set to 1 for 3D psf
%
% OUTPUTS:
%	the updated single molecule parameters
%
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2019.
%	D.Bourgeois, September 2022.
%-

for i=1:n_fluorophores
    n_states=sm_par(i).spectral_data.n_fluorescent_states;

    for j=1:n_states     
        if channel==1
            psf_width=sm_par(i).psf_par_ch1(j).psf_width;
        elseif channel==2
            psf_width=sm_par(i).psf_par_ch2(j).psf_width;
            %Eventually apply defocus in 2D
            if simul_3D==0
                psf_width=psf_width*im_par.two_channel_defocus;
            end
        else
            disp('Error: Unavailable channel number !');
            return;
        end
        
        switch psf_mode
            case 'Simple_Gaussian'
                if simul_3D==0
                    rbox=fix(psf_width/im_par.raster)+1;
                    sig=psf_width/2.35/im_par.raster; % in pixels
                    [xpix,ypix] = meshgrid(-rbox:rbox,-rbox:rbox);
                    psf=exp(-((xpix).*(xpix)+(ypix).*(ypix))/(2*sig^2));
                    %normalize the PSF
                    psf=psf/sum(sum(psf));
                else
                    %calculate a psf kernel over the image z length
                    nslices=im_par.psf_n_zslices;
                    image_height=im_par.raster*im_par.nz; % in [nm]
                    d=im_par.depth_of_focus; % in [nm]
                    zc=im_par.sample_zcenter; % in [nm]
                    if zc==-1 % Case of focal plane at top of sample (for inverted microscope)
                        zc=image_height/2;
                    end
                    c=im_par.nz/2*im_par.raster; % Set focal plane at mid sample height
                    sig=psf_width/2.35/im_par.raster; % in pixels
                    psf(1:nslices)=struct('plane',[],'rbox_x',[],'rbox_y',[]);
                    ABx=im_par.psf_astigmatism_x;
                    ABy=im_par.psf_astigmatism_y;
                    slice_step=image_height/nslices;
                    for k=1:nslices
                        z=k*slice_step+zc; % Should go through the entire image height
                        sig_x=sig*sqrt(max([1e-20,1+((z-c)/d)^2+ABx(1)*((z-c)/d)^3+ABx(2)*((z-c)/d)^4])); % in pixels
                        sig_y=sig*sqrt(max([1e-20,1+((z-c)/d)^2+ABy(1)*((z-c)/d)^3+ABy(2)*((z-c)/d)^4])); % in pixels
                        
                        psf_width_x=2.35*sig_x*im_par.raster;
                        psf_width_y=2.35*sig_y*im_par.raster;
                        
                        rbox_x=fix(psf_width_x/im_par.raster)+1;
                        rbox_y=fix(psf_width_y/im_par.raster)+1;
                        
                        [xpix,ypix]=meshgrid(-rbox_x:rbox_x,-rbox_y:rbox_y);
                        
                        psf_plane=exp(-((xpix).*(xpix))/(2*sig_x^2)-((ypix).*(ypix))/(2*sig_y^2));
                        psf(k).plane=psf_plane/sum(sum(psf_plane));
                        psf(k).rbox_x=rbox_x;
                        psf(k).rbox_y=rbox_y;
                    end
                end
                if simul_3D==0
                    if channel==1
                        sm_par(i).psf_par_ch1(j).psf=psf; % Note that .psf currently not used in get_2Dpeak.m
                        sm_par(i).psf_par_ch1(j).rbox=rbox; 
                    end
                    if channel==2
                        sm_par(i).psf_par_ch2(j).psf=psf; % Note that .psf currently not used in get_2Dpeak.m
                        sm_par(i).psf_par_ch2(j).rbox=rbox; 
                        sm_par(i).psf_par_ch2(j).psf_width=psf_width; % Update psf_width in case of defocus
                    end
                else
                    if channel==1
                        sm_par(i).psf_par_ch1(j).psf=psf; 
                    end
                    if channel==2
                        sm_par(i).psf_par_ch2(j).psf=psf; 
                    end
                end

            otherwise
                disp('PSF Model Not implemented yet')
        end
    end
end
