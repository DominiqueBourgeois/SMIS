function [NbrFit, Fitted_D, Xfit, Yfit] = Fit_SMIS_Lorentzian_Interactive(N, BinEdges, UsedDapp, UseLog, FontSize, x_label, hPlot)

% PURPOSE:
%	Routine to interactively fit diffusion coefficients histograms with Lorentzian models
%
% INPUTS:
%   N: Normalized histogram values
%   BinEdges: Histogram bins edges
%   UsedDapp: The used apparent diffusion coefficients
%   UseLog: Set to 1 if log scale is used
%   FontSize: The font size
%   x_label: The X label
%   hPlot: Handle to the used figure
%
% OUTPUTS:
%   NbrFit: Number of diffusion coefficients fitted
%   Fitted_D: Values of fitted diffusion coefficients
%   Xfit: The fitted X data
%   Yfit: The fitted Y data
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022. Compatibility with SMIS 1.3. Based on
%	software developed by Jean-Bernard Fiche


NbrFit = questdlg('Do you want to fit the distribution with one or two Lorentzian?', 'Fit', '1', '2', '1');
NbrFit = str2double(NbrFit);

BinWidth=BinEdges-circshift(BinEdges,[0,1]);
BinCenters=BinEdges(1:end-1)+BinWidth(2:end)/2;

N=100*N; % Go in percentage
Lorentz = fittype( @(g,x0,A,x) A./(1 + ((x - x0)/g).^2) );

if NbrFit == 1
[A0,Idx] = max(N); % Estimate the starting parameters for the first fit
A0 = double(A0);
X0 = double(BinCenters(Idx));
g0 = double(1);

Xfit = min(BinCenters) : 0.01 : max(BinCenters);
[fitobject,gof1] = fit(double(BinCenters)', double(N)', Lorentz, 'start', [g0,X0,A0]);
Yfit = Lorentz(fitobject.g, fitobject.x0, fitobject.A, Xfit);

figure(hPlot)
ax = gca;
hold off
cla

[N,BinCenters] = hist(UsedDapp, BinCenters);
N_norm = 100*N/sum(N);
b = bar(BinCenters', N_norm);
b(1).FaceColor = [1 0.5 0];

hold on
plot(Xfit, Yfit, '-r', 'LineWidth',1)
ax.FontSize = FontSize;
axis square
box on
xlabel(x_label)
ylabel('Fraction of tracks (%)')


end

if NbrFit == 2
    
    % Fit the two distributions with a Lorentzian function and plot the two
    % fits as well as the sum of the two on the graph. 
    % ------------------------------------------------
    
    warnh = warndlg('Click on the approximate positions of the two maxima.','');
    uiwait(warnh)
    [X,Y] = ginput(2);
    
    X0 = X(1);
    A0 = 100*Y(1);
    X1 = X(2);
    A1 = 100*Y(2);
    
    warnh = warndlg('Click on the approximate positions of the FWHM for the first Lorentzian.','');
    uiwait(warnh)
    [X,~] = ginput(2);
    g0 = abs(X(2)-X(1))/2;
    
    warnh = warndlg('Click on the approximate positions of the FWHM for the second Lorentzian.','');
    uiwait(warnh)
    [X,~] = ginput(2);
    g1 = abs(X(2)-X(1))/2;
    
    Lorentz2 = fittype( @(g1,x01,A1,g2,x02,A2,x) A1./(1 + ((x - x01)/g1).^2) + A2./(1 + ((x - x02)/g2).^2) );
    [fitobject2,gof2] = fit(double(BinCenters)', double(N)', Lorentz2, 'start', [g0,X0,A0,g1,X1,A1]);
    disp(strcat('for the fit, R²=', num2str(100*gof2.rsquare), '%'))
    
    % Find the intersection point between the two fits (estimation)
    % -------------------------------------------------------------
    
    Xfit = min(BinCenters) : 0.01 : max(BinCenters);
    LorentzFit1 = Lorentz(fitobject2.g1, fitobject2.x01, fitobject2.A1, Xfit);
    LorentzFit2 = Lorentz(fitobject2.g2, fitobject2.x02, fitobject2.A2, Xfit);
    
    if fitobject2.x01 < fitobject2.x02
        Idx_Bin = find(Xfit >= fitobject2.x01 & Xfit <= fitobject2.x02);
    else
        Idx_Bin = find(Xfit >= fitobject2.x02 & Xfit <= fitobject2.x01);
    end
    
    Diff = abs(LorentzFit1(Idx_Bin) - LorentzFit2(Idx_Bin));
    [~,Idx] = min(Diff);
    x_Inter = Xfit(Idx_Bin(Idx));
    
    % Analyse separatly the two distributions and return the fraction of
    % events belonging to each of the two distributions. It also splits the
    % trajectories and MSD between the two populations.
    % -------------------------------------------------
       
    if fitobject2.x01 > fitobject2.x02
        Fitted_D = [fitobject2.x01, fitobject2.x02];
        LorentzFit1 = Lorentz(fitobject2.g1, fitobject2.x01, fitobject2.A1, Xfit);
        LorentzFit2 = Lorentz(fitobject2.g2, fitobject2.x02, fitobject2.A2, Xfit);
    else
        Fitted_D = [fitobject2.x02, fitobject2.x01];
        LorentzFit1 = Lorentz(fitobject2.g2, fitobject2.x02, fitobject2.A2, Xfit);
        LorentzFit2 = Lorentz(fitobject2.g1, fitobject2.x01, fitobject2.A1, Xfit);
    end
    Yfit = Lorentz2(fitobject2.g1, fitobject2.x01, fitobject2.A1, fitobject2.g2, fitobject2.x02, fitobject2.A2, Xfit);
    if UseLog==1; Fitted_D=10.^Fitted_D; end
    
    Idx_1 = UsedDapp(:,1) > x_Inter;
    Idx_2 = UsedDapp(:,1) < x_Inter;
    
    UsedDapp_1 = UsedDapp(Idx_1,1);
    UsedDapp_2 = UsedDapp(Idx_2,1);
    F = round(100*size(UsedDapp_1,1)/(size(UsedDapp_1,1)+size(UsedDapp_2,1)));
    disp(strcat('The mobile fraction is equal to f=',num2str(F), '%'));
    
    [N1,BinCenters] = hist(UsedDapp_1(:,1), BinCenters);
    [N2,BinCenters] = hist(UsedDapp_2(:,1), BinCenters);
    N1_norm = 100*N1/sum(N1+N2);
    N2_norm = 100*N2/sum(N1+N2);
    
    % replot the two distributions, showing now the two separate
    % populations
    % -----------
    
    figure(hPlot)
    ax = gca;
    hold off
    cla
    
    b = bar(BinCenters', cat(2, N1_norm', N2_norm'), 'stacked');
    b(1).FaceColor = [1 0.5 0];
    b(2).FaceColor = [0 0.4 1];
    
    hold on
    plot(Xfit, LorentzFit1, '-r', 'LineWidth',1)
    plot(Xfit, LorentzFit2, '-b', 'LineWidth',1)
    plot(Xfit, Yfit, '--k', 'LineWidth',1)
    
    ax.FontSize = FontSize;
    axis square
    box on
    xlabel(x_label)
    ylabel('Fraction of tracks (%)')
    
else
    hold on
    plot(Xfit, Yfit, '--k', 'LineWidth',1)
    disp(strcat('for the fit, R²=', num2str(100*gof1.rsquare), '%'))
    
%     Idx = find(UsedDapp(:,1)>fitobject.x0-3*abs(fitobject.g) & UsedDapp(:,1)<fitobject.x0+3*abs(fitobject.g));
    Fitted_D = fitobject.x0;
    if UseLog==1; Fitted_D=10.^Fitted_D; end
    
    disp(strcat('X0=', num2str(fitobject.x0), '_ s0=', num2str(fitobject.g)));
end

% Display the title on the graph, indicating when there are two populations
% the two diffusion coefficients as well as the fraction of mobile
% particles.
% ----------

    for n = 1 : NbrFit
        
        NewLine = sprintf('D_%d = %.2f µm²/s', n, Fitted_D(1,n));
        if n>1
            Title = strvcat(Title, NewLine);
        else
            Title = NewLine;
        end
        
        if NbrFit>1 && n>1
            NewLine = sprintf('Mobile fraction : %d%%', F);
            Title = strvcat(Title, NewLine);
        end
    end
    
    title(Title);
    
   
