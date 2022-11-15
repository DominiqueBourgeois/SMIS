function sms=trim_sms(sms,im_par)

% PURPOSE:
%	Remove unnecessary fields in sms to save disk space
%
% INPUTS:
%	sms: the single molecules
%   im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the trimmed sms
%
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2021

n_fluorophores=size(sms,2);

%Suppress all the fields that are updated at each frame, so that the final values at the last frame are not interesting
for k=1:n_fluorophores
    sms(k).sm=rmfield(sms(k).sm, 'n_phot_ch1');
    sms(k).sm=rmfield(sms(k).sm, 'n_phot_ch2');
    %     sms(k).sm=rmfield(sms(k).sm, 'n_phot_det_ch1');
    %     sms(k).sm=rmfield(sms(k).sm, 'n_phot_det_ch2');
    sms(k).sm=rmfield(sms(k).sm, 'n_abs');
    sms(k).sm=rmfield(sms(k).sm, 'Ns');
    sms(k).sm=rmfield(sms(k).sm, 'n_em');
    sms(k).sm=rmfield(sms(k).sm, 't_on');
    sms(k).sm=rmfield(sms(k).sm, 'fr_t_on');
    sms(k).sm=rmfield(sms(k).sm, 't_off');
    sms(k).sm=rmfield(sms(k).sm, 'state');
    sms(k).sm=rmfield(sms(k).sm, 'state_id');
    sms(k).sm=rmfield(sms(k).sm, 'state_trace');
    sms(k).sm=rmfield(sms(k).sm, 'fluo_trace');
    sms(k).sm=rmfield(sms(k).sm, 'n_diff_state');
    
    % Clean up the tracks in sptPALM mode
    %We have to take care that the xyz_track coordinates are
    %shifted by one, ie xyz_track(k) correspond to the position BEFORE
    %diffusion is applied at frame k. So that in fact for each frame k
    %the position of the single molecule is xyz_track(k-1)
    if im_par.add_diffusion==1
        for i=1:numel(sms(k).sm)
            if ~isempty(sms(k).sm(i).x_track)
                sms(k).sm(i).x_track(1,:)=[];
                sms(k).sm(i).x_track(:,2)=sms(k).sm(i).x_track(:,2)-1;
                sms(k).sm(i).y_track(1,:)=[];
                sms(k).sm(i).y_track(:,2)=sms(k).sm(i).y_track(:,2)-1;
                if im_par.simul_3D==1
                    sms(k).sm(i).z_track(1,:)=[];
                    sms(k).sm(i).z_track(:,2)=sms(k).sm(i).z_track(:,2)-1;
                end
            end
        end
        sms(k).sm=rmfield(sms(k).sm, 'sub_x');
        sms(k).sm=rmfield(sms(k).sm, 'sub_y');
        if im_par.simul_3D==1
            sms(k).sm=rmfield(sms(k).sm, 'sub_z');
        end    
    end
end