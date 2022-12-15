function [NbrFit, Fitted_D, N, Xfit, Yfit, n_bins, fit_ok] = Fit_D_His_SMIS(Dapp,DiffusionDisplayAndFitMethod)

% PURPOSE:
%	Routine to interactively fit diffusion coefficients histograms with
%	either Gaussian or Lorentzian models, linear or log mode.
%
% INPUTS:
%   Dapp: The apparent diffusion coefficients
%   DiffusionDisplayAndFitMethod: % 1,2: Linear; 3,4: Log10; 1,3: Gaussian fit; 2,4: Lorentzian fit
%
% OUTPUTS:
%   NbrFit: Number of diffusion coefficients fitted
%   Fitted_D: Values of fitted diffusion coefficients
%   N: The normalized histogram values
%   Xfit: The fitted X data
%   Yfit: The fitted Y data
%   n_bins: Used number of bins (50)
%   fit_ok: Yes or No depending if the fit is judged successful
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022. Compatibility with SMIS 1.3. Based on
%	software developed by Jean-Bernard Fiche

%% =======================================================

hPlot=figure; clf; set(gcf, 'color', [1 1 1]); 
% set(gcf, 'Position', get(0, 'Screensize'));
ax = gca;
hold off
cla
FontSize=16;

%Choose between linear and log for diffusion coeffs
LogDapp=log10(Dapp);
if DiffusionDisplayAndFitMethod<3 % Linear 
    UseLog=0;
    UsedDapp = Dapp';
    x_label='Apparent diffusion coefficient (µm²/s)';
elseif DiffusionDisplayAndFitMethod>=3 % Log10
    UseLog=1;
    UsedDapp = LogDapp';
    x_label='Log of apparent diffusion coefficient (µm²/s)';
end


WindowMax=max(UsedDapp(:,1));
WindowMin=min(UsedDapp(:,1));
if isempty(WindowMax)
    WindowMin=0; % Arbitrary number
    WindowMax=1; % Arbitrary number
end
if WindowMax==WindowMin
    WindowMin=WindowMin*0.95;
    WindowMax=WindowMax*1.05;
end

n_bins=50; % Assume 50 different bins
MyStep=(WindowMax-WindowMin)/n_bins;
BinEdges=WindowMin+(0:n_bins)*MyStep; % Produce x_values for histogram
BinEdges(1)=BinEdges(1)-MyStep*1e-6; % Guaranty that no data fall on one of the edge
BinEdges(end)=BinEdges(end)+MyStep*1e-6;


h=histogram(UsedDapp(:,1),'BinEdges',BinEdges,'Normalization','probability');
h.FaceColor = [1 0.5 0];


% [N,Bin] = hist(UsedDapp(:,1), Bin);
% N = 100*N/sum(N);
N=h.Values;
BinEdges=h.BinEdges;

% bar(Bin, N, 'FaceColor', [0 0.4 1])
ax.FontSize = FontSize;
% h.
axis square
box on
xlabel(x_label)
ylabel('Fraction of tracks')


%Choose between Gaussian and Lorentzian for fitting of diffusion coeffs
if DiffusionDisplayAndFitMethod==1 || DiffusionDisplayAndFitMethod==3 % Gaussian fit
    [NbrFit, Fitted_D, Xfit,Yfit] = Fit_SMIS_Gaussian_Interactive(N, BinEdges, UsedDapp, UseLog, FontSize, x_label,hPlot);
elseif DiffusionDisplayAndFitMethod==2 || DiffusionDisplayAndFitMethod==4 % Lorentzian fit
    [NbrFit, Fitted_D, Xfit,Yfit] = Fit_SMIS_Lorentzian_Interactive(N, BinEdges, UsedDapp, UseLog, FontSize, x_label,hPlot);
end
figure(hPlot);
fit_ok=questdlg('Fit OK ?');
% close(hPlot);
