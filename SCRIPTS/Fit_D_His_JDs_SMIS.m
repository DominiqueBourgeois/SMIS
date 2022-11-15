function [NbrFit, Fitted_D, N, Xfit, Yfit, n_bins, fit_ok, Fitted_F] = Fit_D_His_JDs_SMIS(Dapp,tracklength,SX,SY,MyXLim)

% PURPOSE:
%	Routine to interactively fit diffusion coefficients histograms with
%	Model described in Stracy, Mol Cell 81 1499, 2021
%
% INPUTS:
%   Dapp: The apparent diffusion coefficients
%   tracklength= the the length of the tracks that have been used to construct the histogram
%   SX: Starting values for diffusion coefficients
%   SY: Starting values for fractional amplitudes
%   MyXLim: The min,max values for the histogram X-axis
%
% OUTPUTS:
%   NbrFit: Number of diffusion coefficients fitted
%   Fitted_D: Values of fitted diffusion coefficients
%   N: The normalized histogram values
%   Xfit: The fitted X data
%   Yfit: The fitted Y data
%   n_bins: Used number of bins (50)
%   fit_ok: Yes or No depending if the fit is judged successful
%   Fitted_F: Fraction of mobile molecules (1 if 1 diffusion regimes, [immobile,mobile] fraction if 2, [immobile,mobile,most mobile] fractions if 3)
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022. Compatibility with SMIS 1.3. 

%% =======================================================

hPlot=figure; clf; set(gcf, 'color', [1 1 1]); 
ax = gca;
hold off
cla
FontSize=16;

UsedDapp = Dapp';
x_label='Apparent diffusion coefficient (µm²/s)';

if ~exist('MyXLim','var') % If limits not provided
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
else
    WindowMin=MyXLim(1);
    WindowMax=MyXLim(2);
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

ax.FontSize = FontSize;
axis square
box on
xlabel(x_label)
ylabel('Fraction of tracks')

if ~exist('SX','var')
    [NbrFit, Fitted_D, Xfit,Yfit,Fitted_F] = Fit_SMIS_D_From_JDs_Interactive(N, BinEdges, UsedDapp, tracklength, FontSize, x_label,hPlot);
    figure(hPlot);
    fit_ok=questdlg('Fit OK ?');
    
elseif exist('SY','var')
    [NbrFit, Fitted_D, Xfit,Yfit,Fitted_F] = Fit_SMIS_D_From_JDs(N, BinEdges, UsedDapp, tracklength, FontSize, x_label,hPlot,SX,SY);
    fit_ok='Yes'; % Assume the fit will be okay
else
    disp('Insufficient input parameters !')
    NbrFit=[];
    Fitted_D=[];
    Xfit=[];
    Yfit=[];
    Fitted_F=[];
    fit_ok='No';
    return
end

close(hPlot);


