% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&& FITTING OF DATA &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function  [Xfit, Yfit, fitted_values]= fit_histograms(X,Y,...
    fit_type,fit_aniso_type,fit_offset,thermal_relaxation,...
    init, ...
    exp_par,...
    remove_first_point, ...
    fit_data)

global start start_common w_variables w_variables_common nsets datasets_lengths;

nsets=1; % For compatibility with model functions

generate_figure=0;
show_log_scale=0;
%&&&&&&&&&& parameters for fitting &&&&&&&&&&&&&

Xfit=X' ; % my x data
ytofit=Y' ; % my y data

if isfield(init,'a0')
    a0=init.a0;
else
    a0=ytofit(1); % Initial amplitude
end

w_ok=find(ytofit > 0, 1, 'last');
Xfit=Xfit(1:w_ok); % restrict to the non zero values
ytofit=ytofit(1:w_ok);
datasets_lengths=length(Xfit);

if remove_first_point==1     % In this case remove the first point in each sector
    Xfit=Xfit(2:end); % my x data
    ytofit=ytofit(2:end); % my y data
    w_ok=find(ytofit > 0, 1, 'last');
    Xfit=Xfit(1:w_ok); % restrict to the non zero values
    ytofit=ytofit(1:w_ok);
    datasets_lengths=length(Xfit);
end



%Local parameters
%Rates
if thermal_relaxation==0
    epsilon_a=exp_par.epsilon_a; % [mol-1cm-1] exctinction coeff of active molecules (anionic) at used wavelength
    laser_power_density=exp_par.laser_power_density; %[W/cm2] Laser power density at beam center. Get from laser_par.power_density
    laser_wavelength=exp_par.laser_wavelength; %[nm]
    laser_percentage=exp_par.laser_percentage; % [%] Used laser power percentage: Uses output from check_laser_mean_values.m
    laser_i=laser_power_density*laser_wavelength/1.9846e-016; % in ph/cm2/sec
    laser_i=laser_i*laser_percentage/100;
    k_a=laser_i*3.82e-21*epsilon_a; % excitation rate of active molecules
else
    k_a=init.rate1;
    init.phi1=1;
    init.phi2=1;
end


%parameters to fit
%Common parameters
switch fit_type
    case 'monophasic'
        start_common=[init.phi1, 1]';
        fixed_common=[0,  1];
        start = [a0, k_a, init.offset]';
        fixed = [0,     1,      1-fit_offset]; % 1 for fixed parameters
    case 'biphasic'
        if thermal_relaxation==0
            start_common=[init.phi1, init.phi2, 1]';
            fixed_common=[0,  0, 1];
            start = [a0*init.fraction, a0*(1-init.fraction), k_a]';
            fixed = [0,                0,                    1]; % 1 for fixed parameters
        else
            start_common=[init.rate1, init.rate2]'; % phi1:  blinking molecules; phi2: non blinking molecules
            fixed_common=[0,     0];
            start = [a0*init.fraction, a0*(1-init.fraction), init.offset]';
            fixed = [0,                0,                    1-fit_offset]; % 1 for fixed parameters
        end
    case 'stretched'
        start_common=[init.phi1, init.beta]';
        fixed_common=[0,  0];
        start = [a0, k_a, init.offset]';
        fixed = [0,     1,      1-fit_offset]; % 1 for fixed parameters
end



%&&&&&&&&&&&&&&&& END OF INPUT

disp('');
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
disp('Fitting Of Model ...');
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
disp('');

if fit_data==1
    disp('Getting results from fitted values  ...');
end
if ~fit_data
    disp('Getting results from starting values only ...');
end

w_variables_common = find(~fixed_common); n_variables_common=length(w_variables_common);
w_variables = find(~fixed); n_variables=length(w_variables);

if (n_variables+n_variables_common)==0
    disp('All parameters are fixed: cannot fit !');
    return;
end

startvalues_common = start_common(w_variables_common);
startvalues = start(w_variables);
startvalues_all=vertcat(startvalues_common,startvalues);

xdata=Xfit;
ydata=ytofit;

if fit_data==1
    if thermal_relaxation==0
        if strcmp(fit_type,'monophasic')==1 || strcmp(fit_type,'stretched')==1
            if strcmp(fit_aniso_type,'isotropic')
                [fittedvalues_all,r,~,COVB] = nlinfit(xdata, ydata,@fit_exp_function_multi,startvalues_all);
            end
            if strcmp(fit_aniso_type,'aniso')
                [fittedvalues_all,r,~,COVB] = nlinfit(xdata, ydata,@fit_exp_function_multi_aniso,startvalues_all);
            end
            if strcmp(fit_aniso_type,'aniso_linear')
                [fittedvalues_all,r,~,COVB] = nlinfit(xdata, ydata,@fit_exp_function_multi_aniso_linear,startvalues_all);
            end
        end
        
        if strcmp(fit_type,'biphasic')==1
            if strcmp(fit_aniso_type,'isotropic')
                [fittedvalues_all,r,~,COVB] = nlinfit(xdata, ydata,@fit_biexp_function_multi,startvalues_all);
            end
            if strcmp(fit_aniso_type,'aniso')
                [fittedvalues_all,r,~,COVB] = nlinfit(xdata, ydata,@fit_biexp_function_multi_aniso,startvalues_all);
            end
            if strcmp(fit_aniso_type,'aniso_linear')
                [fittedvalues_all,r,~,COVB] = nlinfit(xdata, ydata,@fit_biexp_function_multi_aniso_linear,startvalues_all);
            end
        end
    end
    
    if thermal_relaxation==1
        if strcmp(fit_type,'monophasic')==1 || strcmp(fit_type,'stretched')==1
            [fittedvalues_all,r,~,COVB] = nlinfit(xdata, ydata,@fit_exp_function_multi,startvalues_all);
        end
        if strcmp(fit_type,'biphasic')==1
            [fittedvalues_all,r,~,COVB] = nlinfit(xdata, ydata,@fit_2_rates_multi,startvalues_all);
        end
    end
    
    [ci, ~] = nlparci_dom(fittedvalues_all,r,'covar',COVB); %calculate error estimates on parameters
    ci_err=(ci(:,2)-ci(:,1))/2; % 95% confidence Standard error is in 'se' variable
    
    
    result_common = start_common;
    fit_error_common=zeros(size(result_common));
    result = start;
    fit_error=zeros(size(result));
    if ~isempty(w_variables_common)
        result_common(w_variables_common) = fittedvalues_all(1:length(w_variables_common));
        fit_error_common(w_variables_common) = ci_err(1:length(w_variables_common));
    end
    if ~isempty(w_variables)
        result(w_variables) = fittedvalues_all(length(w_variables_common)+1:end);
        fit_error(w_variables) = ci_err(length(w_variables_common)+1:end);
    end
else
    result_common = start_common;
    fit_error_common=zeros(size(result_common));
    result = start;
    fit_error=zeros(size(result));
end


%Get back the results
result=[result',result_common'];
fit_error=[fit_error',fit_error_common'];

if thermal_relaxation==0
    if strcmp(fit_type,'monophasic')==1 || strcmp(fit_type,'stretched')==1
        fitted_values.fitted_a0=result(1);
        fitted_values.fitted_rate=result(2)*result(4); % Observed rate
        fitted_values.fitted_offset=result(3);
        fitted_values.fitted_yield=result(4);
        fitted_values.fitted_beta=result(5);
        
        fitted_values.fitted_a0_error=fit_error(1);
        fitted_values.fitted_rate_error=result(2)*fit_error(4);
        fitted_values.fitted_offset_error=fit_error(3);
        fitted_values.fitted_yield_error=fit_error(4);
        fitted_values.fitted_beta_error=fit_error(5);
        
        if strcmp(fit_aniso_type,'isotropic')
            Yfit = get_exp_model(result, Xfit); % Get yfit
        end
        if strcmp(fit_aniso_type,'aniso')
            Yfit = get_exp_aniso_model(result, Xfit); % Get yfit
        end
        if strcmp(fit_aniso_type,'aniso_linear')
            Yfit = get_exp_aniso_linear_model(result, Xfit); % Get yfit
        end
    end
    
    if strcmp(fit_type,'biphasic')==1
        
        fitted_values.fitted_a0_1=result(1);
        fitted_values.fitted_a0_2=result(2);
        fitted_values.fitted_rate_1=result(3)*result(4); % Observed rate 1
        fitted_values.fitted_rate_2=result(3)*result(5); % Observed rate 2
        fitted_values.fitted_yield_1=result(4);
        fitted_values.fitted_yield_2=result(5);
        fitted_values.fitted_beta=result(6);       
        
        fitted_values.fitted_a0_1_error=fit_error(1);
        fitted_values.fitted_a0_2_error=fit_error(2);
        fitted_values.fitted_rate_1_error=result(3)*fit_error(4);
        fitted_values.fitted_rate_2_error=result(3)*fit_error(5);
        fitted_values.fitted_yield_1_error=fit_error(4);
        fitted_values.fitted_yield_2_error=fit_error(5);
        fitted_values.fitted_beta_error=fit_error(6);
        
        if strcmp(fit_aniso_type,'isotropic')
            Yfit = get_biexp_model(result, Xfit); % Get yfit
        end
        if strcmp(fit_aniso_type,'aniso')
            Yfit = get_biexp_aniso_model(result, Xfit); % Get yfit
        end
        if strcmp(fit_aniso_type,'aniso_linear')
            Yfit = get_biexp_aniso_linear_model(result, Xfit); % Get yfit
        end
    end
end




if thermal_relaxation==1
    if strcmp(fit_type,'monophasic')==1 || strcmp(fit_type,'stretched')==1
        
        fitted_values.fitted_a0=result(1);
        fitted_values.fitted_rate=result(2)*result(4); % Observed rate
        fitted_values.fitted_offset=result(3);
        fitted_values.fitted_beta=result(5);
        
        fitted_values.fitted_a0_error=fit_error(1);
        fitted_values.fitted_rate_error=result(2)*fit_error(4);
        fitted_values.fitted_offset_error=fit_error(3);
        fitted_values.fitted_beta_error=fit_error(5);     
        
        Yfit = get_exp_model(result, Xfit); % Get yfit
    end
    if strcmp(fit_type,'biphasic')==1
        fitted_values.fitted_a0_1=result(1);
        fitted_values.fitted_a0_2=result(2);
        fitted_values.fitted_offset=result(3);
        fitted_values.fitted_rate_1=result(4); % Observed rate 1
        fitted_values.fitted_rate_2=result(5); % Observed rate 1
                
        fitted_values.fitted_a0_1_error=fit_error(1);
        fitted_values.fitted_a0_2_error=fit_error(2);
        fitted_values.fitted_offset_error=fit_error(3);
        fitted_values.fitted_rate_1_error=result(3)*fit_error(4);
        fitted_values.fitted_rate_2_error=result(3)*fit_error(5);       
        
        Yfit = get_2_rates_model(result, Xfit); % Get yfit
    end
end


if generate_figure==1
    figure('Color',[1 1 1]);
    fig_fontsize=16;
    hold on
    if remove_first_point==1 % make sure datasets_lengths are correct
        datasets_lengths=datasets_lengths+1;
    end
    
    bar(Xfit, ytofit, 0.75, 'c');
    if remove_first_point==0
        line(Xfit,  Yfit,'Color','Blue','LineWidth',3);
    else
        line(Xfit,  Yfit,'Color','Green','LineWidth',3);
    end
    
       
    if show_log_scale==1
        set(gca,'YScale','log')
        tmp_ylim=ylim;
        ylim([0.1, tmp_ylim(2)]);
    end
    xtitle='Time [s]';
    ytitle='# Mol';
    xlabel(xtitle,'fontsize',fig_fontsize,'fontweight','b')
    ylabel(ytitle,'fontsize',fig_fontsize,'fontweight','b')
end

%Final comparison of initial data series and fitted data series

% disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
% disp('&&&&& STARTING VALUES &&&&&&');
% disp(start);
% disp('&&&&& FITTED VALUES &&&&&&');
% disp(result);
% disp('&&&&&  AMPLITUDES &&&&&&');
% disp(fitted_a0);
% disp('&&&&& ERROR (95% CONFIDENCE) &&&&&&');
% disp(fitted_a0_error);
% disp('&&&&&  RATES [s-1] &&&&&&');
% disp(fitted_rates*fitted_yield);
% disp('&&&&& ERROR (95% CONFIDENCE) &&&&&&');
% disp(fitted_rates_error*fitted_yield_error);
% disp('&&&&&  YIELD &&&&&&');
% disp(fitted_yield);
% disp('&&&&& ERROR (95% CONFIDENCE) &&&&&&');
% disp(fitted_yield_error);
% if fit_beta==1
%     disp('&&&&&  BETA (stretched exponential) &&&&&&');
%     disp(fitted_beta);
%     disp('&&&&& ERROR (95% CONFIDENCE) &&&&&&');
%     disp(fitted_beta_error);
% end
% if fit_offset==1
%     disp('&&&&&  OFFSET &&&&&&');
%     disp(fitted_offset);
%     disp('&&&&& ERROR (95% CONFIDENCE) &&&&&&');
%     disp(fitted_offset_error);
% end
% 
disp('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');

fitted_values.result=result; 
fitted_values.fit_error=fit_error;

