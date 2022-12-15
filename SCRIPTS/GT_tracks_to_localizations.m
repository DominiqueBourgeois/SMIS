function [l] = GT_tracks_to_localizations(t,raster)

% Converts a track structure towards a localisation structure
% for each localization: x, y, f: (frame #), tn: track number, tn2:
% matching track number (from exp or GT), m: 1 if matched 

n_t=length(t);

k=1;
% Data x,y in raster, convert to nm
n_l=sum([t.tracklength]);
l(1:n_l)=struct('x',0, 'y',0,'f',0,'tn',0,'tn2',[],'m',0);

% if mode==0 % GT data
    for i=1:n_t
        newVals = num2cell([t(i).x]*raster); [l(k:k-1+t(i).tracklength).x] = newVals{:};
        newVals = num2cell([t(i).y]*raster); [l(k:k-1+t(i).tracklength).y] = newVals{:};
        newVals = num2cell([t(i).frames]); [l(k:k-1+t(i).tracklength).f] = newVals{:};
        newVals = num2cell(i); [l(k:k-1+t(i).tracklength).tn] = deal(newVals{:});
        k=k+t(i).tracklength;
    end
% end

% if mode==1 % SWIFT data 
%     n_l=sum([t.tracklength]);
%     l(1:n_l)=struct('x',0, 'y',0,'f',0,'tn',0,'tn2',[],'m',[]);
% 
%     for i=1:n_t
%         newVals = num2cell([t(i).x]); [l(k:k-1+t(i).tracklength).x] = newVals{:};
%         newVals = num2cell([t(i).y]); [l(k:k-1+t(i).tracklength).y] = newVals{:};
%         newVals = num2cell([t(i).frames]); [l(k:k-1+t(i).tracklength).f] = newVals{:};
%         newVals = num2cell(i); [l(k:k-1+t(i).tracklength).tn] = deal(newVals{:});
%         k=k+t(i).tracklength;
%     end
% end


end

