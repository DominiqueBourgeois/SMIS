function im=position_peak(spot,xy,im)
% NAME:
%	POSITION_PEAK
%
% PURPOSE:
%	Place a single molecule spot into a large image
%
% CATEGORY:
%	Signal, image processing.
%
% CALLING SEQUENCE:
%	position_peak(spot,xy,im);
%
% INPUTS:
%	spot: the single molecules spot (small array)
%	xy:	the central position of the spot in image
%   im: the image to be modified (large array)
%
% OUTPUTS:
%	im: the modified image where the pixel intensities of spot have been
%	added.
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
%	D.Bourgeois, September 2020. Do nothing if spot is empty. Separate x
%	and y

if ~isempty(spot)
    x=xy(1); % spot centered on pixel # containing the molecule
    y=xy(2);
    ss=size(spot);
    rbox=fix(ss/2);
    xs=size(im,1);
    ys=size(im,2);
    
    xr=(x-rbox(1):x+rbox(1));
    yr=(y-rbox(2):y+rbox(2));
    
    if min(xr)>0 && max(xr)<=xs && min(yr)>0 && max(yr)<=ys
        im(xr,yr)=im(xr,yr)+spot;
    else % Treat case where spot not within FOV
        if xr(1)<=xs && yr(1)<=ys && xr(end)>=1 && yr(end)>=1 % See if part of spot within image
            
            wxrs=find(xr>0);
            xstart=wxrs(1);
            
            wxre=find(xr>=xs);
            if isempty(wxre)
                xend=ss(1);
            else
                xend=wxre(1);
            end
            
            wyrs=find(yr>0);
            ystart=wyrs(1);
            
            wyre=find(yr>=ys);
            if isempty(wyre)
                yend=ss(2);
            else
                yend=wyre(1);
            end
            
            reduced_spot=spot((xstart:xend),(ystart:yend));
            xss=xstart-1;
            xes=ss(1)-xend;
            yss=ystart-1;
            yes=ss(2)-yend;
            
            xrr=(x-rbox(1)+xss:x+rbox(1)-xes);
            yrr=(y-rbox(2)+yss:y+rbox(2)-yes);
            im(xrr,yrr)=im(xrr,yrr)+reduced_spot;
        end
    end
end


