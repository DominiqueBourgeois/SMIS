function sm_par=get_psf_width(n_fluorophores,sm_par,im_par)

% PURPOSE:
%	Get width of PSF
%
% INPUTS:
%   n_fluorophores: the # of fluorophores
%	sm_par: the sm parameters
%
% OUTPUTS:
%	sm_par: the updated sm parameters for psf_width
%
% MODIFICATION HISTORY:
%	D.Bourgeois, October 2019: version > simulate_palm_vsn15
% 	D.Bourgeois, September 2022: change filters definitions

obj_na=im_par.obj.na;

for i=1:n_fluorophores
    n_states=sm_par(i).spectral_data.n_fluorescent_states;
    psf_par_ch1(1:n_states)=struct('psf_width',[],'psf',[], 'rbox',[]);
    if im_par.two_channel==1
        psf_par_ch2(1:n_states)=struct('psf_width',[],'psf',[], 'rbox',[]);
    end

    for j=1:n_states
        em_spectrum=sm_par(i).spectral_data.em_spectra(j).s;
        s=em_spectrum(:,2); % The spectral values
        l=em_spectrum(:,1); % The wavelength values

        %Channel 1
        channel=1;
        f=sm_par(i).filter_profiles(j).ch1; % Filter profile for ch1

        if im_par.two_channel==1 && im_par.filters.add_dichroic==1
            df=sm_par(i).filter_profiles(j).df; % The dichroic
            fs=s.*f.*(1-df);
        else 
            fs=s.*f;
        end

        if sum(fs)>0 % Check that filtered spectrum is nonzero
            lambda=sum(fs.*l)/sum(fs);
            psf_par_ch1(j).psf_width = 1.22*lambda/(2*obj_na);
        else
            psf_par_ch1(j).psf_width=0 ;
            disp(['**** Fluorophore: ',sm_par.fluorophore_name,' has no emission in channel ',num2str(channel),' for fluorescent state ', num2str(j),' ; setting arbitrary PSF width to 0 nm !']);
        end
        disp([sm_par(i).fluorophore_name,': PSF FWHM for fluorescent state ', num2str(j),', Channel ',num2str(channel),'  [nm]: ',num2str(psf_par_ch1(j).psf_width)]);

        %Channel 2
        if im_par.two_channel==1
            channel=2;
            f=sm_par(i).filter_profiles(j).ch2; % Filter profile for ch2

            if im_par.filters.add_dichroic==1
                df=sm_par(i).filter_profiles(j).df; % The dichroic
                fs=s.*f.*df;
            else
                fs=s.*f;
            end

            if sum(fs)>0 % Check that filtered spectrum is nonzero
                lambda=sum(fs.*l)/sum(fs);
                psf_par_ch2(j).psf_width = 1.22*lambda/(2*obj_na);
            else
                psf_par_ch2(j).psf_width=0 ;
                disp(['**** Fluorophore: ',sm_par.fluorophore_name,' has no emission in channel ',num2str(channel),' for fluorescent state ', num2str(j),' ; setting arbitrary PSF width to 0 nm !']);
            end
            disp([sm_par(i).fluorophore_name,': PSF FWHM for fluorescent state ', num2str(j),', Channel ',num2str(channel),'  [nm]: ',num2str(psf_par_ch2(j).psf_width)]);
        end

    end

    sm_par(i).psf_par_ch1=psf_par_ch1;
    if im_par.two_channel==1
        sm_par(i).psf_par_ch2=psf_par_ch2;
    end

end