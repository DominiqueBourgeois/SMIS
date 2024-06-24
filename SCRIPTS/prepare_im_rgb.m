function true_im_rgb=prepare_im_rgb(a_h_all, im_par)

% PURPOSE:
%	Prepare a rgb image with all true positions of all sms
%
% INPUTS:
%	a_h_all: the images onto place the SMs
%	im_par: general parameters on the imaging experiment
%
% OUTPUTS:
%   true_im_rgb = the RGB image, where the SMs have been positionned
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15

n_fluorophores=size(a_h_all, 1);
% if n_fluorophores > 6 
%     error('prepare_im_rgb.m: Not enough colors defined for more than 6 dyes ...')
% end

RGB=zeros(3,max([6,n_fluorophores]));
RGB(:,1)=[0,1,0]; % Green
RGB(:,2)=[1,0,0]; % Red
RGB(:,3)=[0,0,1]; % Blue
RGB(:,4)=[1,1,0]; % Yellow
RGB(:,5)=[1,0,1]; % Magenta
RGB(:,6)=[0,1,1]; % Cyan
if n_fluorophores>6
    RGB(:,7:end)=rand(3,n_fluorophores-6); % Random color
end

true_im_rgb=zeros(im_par.binning*im_par.n,im_par.binning*im_par.m,3);
for i=1:n_fluorophores
    true_im_rgb(:,:,1)=true_im_rgb(:,:,1)+squeeze(RGB(1,i)*a_h_all(i,:,:)); % set to red RGB=[0,1,0]
    true_im_rgb(:,:,2)=true_im_rgb(:,:,2)+squeeze(RGB(2,i)*a_h_all(i,:,:)); % set to green RGB=[0,1,0]
    true_im_rgb(:,:,3)=true_im_rgb(:,:,3)+squeeze(RGB(3,i)*a_h_all(i,:,:)); % set to blue RGB=[0,1,0]   
end

% Normalize image
true_im_rgb(:,:,1)=true_im_rgb(:,:,1)./max(max(true_im_rgb(:,:,1)));
true_im_rgb(:,:,2)=true_im_rgb(:,:,2)./max(max(true_im_rgb(:,:,2)));
true_im_rgb(:,:,3)=true_im_rgb(:,:,3)./max(max(true_im_rgb(:,:,3)));



end