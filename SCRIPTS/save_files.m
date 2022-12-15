function save_files(simulation_files,im_par,sm_par,sms,lasers)

disp('Saving data ...');
%Write out the true high-resolution image
disp('Saving true image ...');
if im_par.simul_3D==0
    imwrite_dom(im_par.true_im_rgb,simulation_files.im_true_out,'tif');
elseif im_par.simul_3D==1
    true_3D_kernel=im_par.true_3D_kernel;
    save(simulation_files.im_true_out_3D,'true_3D_kernel');
end
%Write out the diffraction limited (sum) images
if max(max(im_par.emccd_im_ch1_dl))>=2^16 % Normalize if necessary to 16 bits
    im_par.emccd_im_ch1_dl=im_par.emccd_im_ch1_dl*(2^16-1)/max(max(im_par.emccd_im_ch1_dl));
end
if im_par.two_channel==1
    if max(max(im_par.emccd_im_ch2_dl))>=2^16 % Normalize if necessary to 16 bits
        im_par.emccd_im_ch2_dl=im_par.emccd_im_ch2_dl*(2^16-1)/max(max(im_par.emccd_im_ch2_dl));
    end
end
disp('Saving diffraction limited image (channel 1) ...');
imwrite_dom(uint16(im_par.emccd_im_ch1_dl),simulation_files.im_ch1_dl_out,'tif');
if im_par.two_channel==1
    disp('Saving diffraction limited image (channel 2) ...');
    imwrite_dom(uint16(im_par.emccd_im_ch2_dl),simulation_files.im_ch2_dl_out,'tif');
end
disp('Saving laser profile ...');
fields = fieldnames(simulation_files);
n_lasers=size(lasers,2);
for k=1:n_lasers
    k2=size(fields,1)-(n_lasers-k);
    imwrite_dom(uint16(lasers(k).beam_profile*65535/max(max(lasers(k).beam_profile))),simulation_files.(fields{k2}),'tif');
end
disp('Saving single molecule data ...');
for i=1:size(sm_par,2)
    sm_par(i).w_patterns=[]; % Do not save this to save time & space
end
im_par.TiffInfo=[]; % We don't want this in the saved data (generate a warning message upon loading)
save(simulation_files.data_out,'im_par','sm_par','sms','lasers'); %Save data
disp('Good, your SMIS Data are saved !');
end




