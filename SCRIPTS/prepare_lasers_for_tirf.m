function lasers=prepare_lasers_for_tirf(n_lasers, lasers, im_par)

%
% PURPOSE:
%   Prepare lasers for TIRF or HILO mode
%
% INPUTS:
%   n_lasers: # of lasers
%   lasers: lasers to analyse
%	im_par: the imaging parameters
%
% OUTPUTS:
%   lasers: updated lasers
%
% MODIFICATION HISTORY:
%	D.Bourgeois, November 2020.


if im_par.simul_3D==1 && any([lasers.tirf]==1)
    disp(['TIRF Critical angle is: ',num2str(im_par.obj.critical_angle*180/pi),'°']);
    beta=im_par.obj.sample_indice/im_par.obj.immersion_indice;
    for i=1:n_lasers % set tirf angle for automatized choice
        if lasers(i).tirf_angle==-1
            lasers(i).tirf_angle=1.05*180/pi*im_par.obj.critical_angle; % Use 5% more than theoritical critical angle
        end
        theta=lasers(i).tirf_angle*pi/180; % theta angle in rad
        
        if theta>=im_par.obj.critical_angle   % TIRF mode

            % Intensity for p-polarized part Axelrod Meth Enz 2003 eq 10
            Ip=4*cos(theta)^2*(2*sin(theta)^2-beta^2)/(beta^4*cos(theta)^2+sin(theta)^2-beta^2);
            % Intensity for s-polarized part Axelrod Meth Enz 2003 eq 11
            Is=4*cos(theta)^2/(1-beta^2);
            %For a circularly polarized beam the contributions of Ip and Is are
            %equal
            lasers(i).tirf_amplification=(Ip+Is)/2;
            
            %Set the characteristic depth for TIRF penetration
            lasers(i).d=lasers(i).wavelength/(4*pi)/sqrt(im_par.obj.immersion_indice^2*(sin(theta))^2-im_par.obj.sample_indice^2); % in [nm]
        else % HILO mode
            lasers(i).tirf_amplification=1; % No amplification in this case
        end
    end
end





