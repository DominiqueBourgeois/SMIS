function [fitted_values,fit_ok, Xfit, Yfit]=Fit_MyData(h, fit_type, remove_first_point, fit_data)

% PURPOSE:
%	Routine to fit data with various flavors of exponential decay's
%
% INPUTS:
%   h: Histogram data
%   fit_type: The model function to use (exponential, biexponential, stretched exponential)
%   remove_first_point: Set to 1 to remove the first experimental point from the fit
%   fit_data: Set to 1 (or 0) to actually perform (or not) the fit
%
% OUTPUTS:
%   fitted_values: The fitted model parameters
%   fit_ok: 1 if fit was successful
%   Xfit: The fitted X data
%   Yfit: The fitted Y data
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022. Compatibility with SMIS 1.3. 


MyBinWidth=h.BinEdges-circshift(h.BinEdges,[0,1]);
MyBinWidth=MyBinWidth(2:end);

X=h.BinEdges(1:end-1)+MyBinWidth/2.;
Y=h.Values;

fit_aniso_type=0;
fit_offset=0;
thermal_relaxation=1;

init.a0=Y(1);
init.phi1=nan;
init.phi2=nan;
init.beta=1;
%estimate initial rate
t=abs(Y/Y(1)-1/2);
w=find(t==min(t),1);
init.rate1=log(2)/(X(w)-X(1));
init.rate2=init.rate1/2; % Usually works fine
init.offset=0;
init.fraction=0.75; % Initial fraction of population 1, usually works fine
exp_par=[];

fit_ok=1;

lastwarn('') % Clear last warning message
try
    [Xfit, Yfit, fitted_values]=...
        fit_histograms(...
        X,Y,...
        fit_type,fit_aniso_type,fit_offset,thermal_relaxation,...
        init, ...
        exp_par,...
        remove_first_point, ...
        fit_data);
catch
    disp('Could not fit the data ! Consider changing fitting method !');
    fit_ok=0;
end

[warnMsg, ~] = lastwarn;
if ~isempty(warnMsg)
    disp('Could not properly fit the data ! Consider changing fitting method !');
    fit_ok=0.5;
end

% deep_blue=[40 35 225]/255;
% deep_yellow=[220 240 40]/255;
% yellow=[255 255 0]/255;
% purple=[128 0 128]/255;
% magenta=[255 0 255]/255;
% if fit_ok>0 % Plot the data if fit is ok
%     if remove_first_point==0
%         line(ax,Xfit,  Yfit,'Color','Blue','LineWidth',3);
%     else
%         line(ax,Xfit,  Yfit,'Color','Green','LineWidth',3);
%     end
% else
%     fitted_values=nan;
% end

if fit_ok==false
    fitted_values=nan;
end
    
