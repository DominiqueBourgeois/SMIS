function R=fit_SMIS_MSDMeans(Lags,Msds,Msds_Std, par)

% PURPOSE:
%	Fit the MSD's with the model: MSDs=2d*D*t^alpha+offset, where d is the
%   dimensionality (2 for 2D processing, 3 3D processing), alpha is a
%   coefficient (=1 for free diffusion, <1 for constrainted diffusion, > 1
%   for super diffusive motion), and offset is the offset.
%
% INPUTS:
%   Lags: The lags [s]
%   Msds: The corresponding mean MSDs [um²]
%   Msds_Std: The corresponding standard deviations of the MSDs [um²]
%   par: Input parameters
%       3D_Mode: set to 1 for 3D fit
%
% OUTPUTS:
%   R: Structure containing fitted D, alpha, offset as well as standard errors
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022. Compatibility with SMIS 1.3

% First do a simple linear effect to estimate the starting values
xdata=Lags;
ydata=Msds;

Weight=1./Msds_Std';
Weight(isinf(Weight))=0.1*min(Weight(~isinf(Weight))); % trick to avoid +inf Weights

[fitobject,~] = fit(xdata', ydata', 'Poly1', 'Weights', Weight, 'Upper', [Inf, min(ydata)]);
Dstart=fitobject.p1/4; % Return the coefficient of diffusion

% Define lower and upper bounds
lb=[0 0 -Inf];
ub=[Inf Inf Inf];

%Eventually fix Alpha
if par.C(1)==1
    lb(2)=1;
    ub(2)=1;
end

%Eventually fix Offset
if par.C(2)==1
    lb(3)=0;
    ub(3)=0;
end

% Define fit options
fit_options = optimoptions('lsqcurvefit');
% fit_options.PlotFcn='optimplotresnorm';

% Define fitting function Y=2*d*D*t^alpha+offset
fh = @(x,xdata)(2*par.d*x(1)*xdata.^x(2)+x(3));

try
    startvalues=[Dstart,1,0]; % Start with alpha=1 and offset=0;
    [fittedvalues, ~, r, ~, ~, ~, J] = lsqcurvefit(fh,startvalues, xdata, ydata, lb, ub, fit_options);
catch
    warndlg('Fitting failed ! ');
    R=[];
    return
end

[ci, ~] = nlparci_SMIS(fittedvalues,r,'jacobian',J); %calculate error estimates on parameters
ci_err=(ci(:,2)-ci(:,1))/2; % 95% confidence Standard error is in 'se' variable

R.D=fittedvalues(1);
R.alpha=fittedvalues(2);
R.offset=fittedvalues(3);
R.D_error=ci_err(1);
R.alpha_error=ci_err(2);
R.offset_error=ci_err(3);



