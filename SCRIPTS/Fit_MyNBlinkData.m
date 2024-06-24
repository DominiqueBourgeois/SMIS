% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&& FITTING OF DATA &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function [fitted_eta, fitted_eta_error, fit_ok, xtofit, yfit]=Fit_MyNBlinkData(...
    h,initial_eta, remove_first_point, fit_data)

global start w_variables;

%&&&&&&&&&& parameters for fitting &&&&&&&&&&&&&
% eta=kd/(kd+kb) kd=rate for blinking off; kb=bleaching rate

xtofit=h.BinEdges(1:end-1)+h.BinWidth/2; % my x data
ytofit=h.Values; % my y data
initial_a0=ytofit(1);

%fit the data ? (1: yes% 0: no, just apply the starting model)
fit = fit_data;

%parameters to fit
%Common parameters

start=[initial_a0, initial_eta]';
fixed=[0,  0];

%&&&&&&&&&&&&&&&& END OF INPUT

disp('');
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
disp('Fitting Of Model ...');
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
disp('');

if fit==1
    disp('Getting results from fitted values  ...');
end
if ~fit
    disp('Getting results from starting values only ...');
end

w_variables = find(~fixed); n_variables=length(w_variables);

if n_variables==0
    disp('All parameters are fixed: cannot fit !');
    return;
end

startvalues = start(w_variables);
startvalues_all=startvalues;


%Possibly do the fit on the blinking molecules only, by removing the first
%data point (corresponding to non blinking molecules): do this if 2
%populations of molecules are present: the blinking ones and the
%non-blinking ones
if remove_first_point==1 && length(xtofit)>1
    xdata=xtofit(2:end);
    ydata=ytofit(2:end);
else
    xdata=xtofit;
    ydata=ytofit;
end

fit_ok=1;

if fit==1
    
    lastwarn('') % Clear last warning message
    try
        [fittedvalues_all,r,~,COVB] = nlinfit(xdata, ydata,@fit_N_blink_hist,startvalues_all);
    catch
        disp('Could not fit the N-blink data ! Consider changing initial p-value ? ');
        fit_ok=0;
    end
    
    [warnMsg, ~] = lastwarn;
    if ~isempty(warnMsg)
        disp('Could not properly fit the N_blink data !');
        fit_ok=0.5;
    end
    
    if fit_ok>0
        
        [ci, ~] = nlparci_dom(fittedvalues_all,r,'covar',COVB); %calculate error estimates on parameters
        ci_err=(ci(:,2)-ci(:,1))/2; % 95% confidence Standard error is in 'se' variable
        
        result = start;
        fit_error=zeros(size(result));
        
        if ~isempty(w_variables)
            result(w_variables) = fittedvalues_all(1:length(w_variables));
            fit_error(w_variables) = ci_err(1:length(w_variables));
        end
    else
        result = nan;
        fit_error=nan(size(result));
    end
else
    result = start;
    fit_error=zeros(size(result));
end

if fit_ok>0
    fitted_a0=result(1);
    fitted_a0_error=fit_error(1);
    fitted_eta=result(2);
    fitted_eta_error=fit_error(2);
    
    yfit = get_N_blink_model(result, xtofit); % Get yfit
 else
    fitted_eta=nan;
    fitted_eta_error=nan;
    fitted_a0=nan;
    fitted_a0_error=nan;
end
    
    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
    disp('&&&&&  AMPLITUDES &&&&&&');
    disp(fitted_a0);
    disp('&&&&& ERROR (95% CONFIDENCE) &&&&&&');
    disp(fitted_a0_error);
    disp('&&&&&  p-value &&&&&&');
    disp(1-fitted_eta);
    disp('&&&&& ERROR (95% CONFIDENCE) &&&&&&');
    disp(fitted_eta_error);
    
    disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
end

