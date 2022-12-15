function [NbrFit, Fitted_D, Xfit, Yfit] = Fit_SMIS_Gaussian_Interactive(N, BinEdges, UsedDapp, UseLog, FontSize, x_label, hPlot)

% PURPOSE:
%	Routine to interactively fit diffusion coefficients histograms with Gaussian models
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




NbrFit = questdlg('Do you want to fit the distribution with one, two or three Gaussian ?', 'Fit', '1', '2','3','1');
NbrFit = str2double(NbrFit);

BinWidth=BinEdges-circshift(BinEdges,[0,1]);
BinCenters=BinEdges(1:end-1)+BinWidth(2:end)/2;

N=100*N; % Go in percentage
Gauss = fittype( @(sig,x0,A,x) A.*(exp(-((x - x0).^2)./(2*sig^2))));

if NbrFit == 1
    [A0,Idx] = max(N); % Estimate the starting parameters for the first fit
    A0 = double(A0);
    X0 = double(BinCenters(Idx));
    sig0 = double(0.1);
    Xfit = min(BinCenters) : 0.01 : max(BinCenters);
    [fitobject,~] = fit(double(BinCenters)', double(N)', Gauss, 'start', [sig0,X0,A0]);
    Yfit = Gauss(fitobject.sig, fitobject.x0, fitobject.A, Xfit);
    Fitted_D = fitobject.x0;
    if UseLog==1; Fitted_D=10.^Fitted_D; end
    
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
    
    % Fit the two distributions with a Gaussian function and plot the two
    % fits as well as the sum of the two on the graph. 
    % ------------------------------------------------
    
    warnh = warndlg('Click on the approximate positions of the two maxima.','');
    uiwait(warnh)
    [X,Y] = ginput(2);
    
    X0 = X(1);
    A0 = 100*Y(1);
    X1 = X(2);
    A1 = 100*Y(2);
    
    warnh = warndlg('Click on the approximate positions of the FWHM for the first Gaussian.','');
    uiwait(warnh)
    [X,~] = ginput(2);
    sig0 = abs(X(2)-X(1))/2.35;
    
    warnh = warndlg('Click on the approximate positions of the FWHM for the second Gaussian.','');
    uiwait(warnh)
    [X,~] = ginput(2);
    sig1 = abs(X(2)-X(1))/2.35;
    
    Gauss2 = fittype( @(sig1,x01,A1,sig2,x02,A2,x) A1.*(exp(-((x - x01).^2)./(2*sig1^2))) + A2.*(exp(-((x - x02).^2)./(2*sig2^2))) );
    [fitobject2,gof2] = fit(double(BinCenters)', double(N)', Gauss2, 'start', [sig0,X0,A0,sig1,X1,A1]);
    disp(strcat('for the fit, R²=', num2str(100*gof2.rsquare), '%'))
    
    % Find the intersection point between the two fits (estimation)
    % -------------------------------------------------------------
    
    Xfit = min(BinCenters) : 0.01 : max(BinCenters);
    GaussFit1 = Gauss(fitobject2.sig1, fitobject2.x01, fitobject2.A1, Xfit);
    GaussFit2 = Gauss(fitobject2.sig2, fitobject2.x02, fitobject2.A2, Xfit);
    
    if fitobject2.x01 < fitobject2.x02
        Idx_Bin = find(Xfit >= fitobject2.x01 & Xfit <= fitobject2.x02);
    else
        Idx_Bin = find(Xfit >= fitobject2.x02 & Xfit <= fitobject2.x01);
    end
    
    Diff = abs(GaussFit1(Idx_Bin) - GaussFit2(Idx_Bin));
    [~,Idx] = min(Diff);
    x_Inter = Xfit(Idx_Bin(Idx));
    
    % Analyse separatly the two distributions and return the fraction of
    % events belonging to each of the two distributions. It also splits the
    % trajectories and MSD between the two populations.
    % -------------------------------------------------
       
    if fitobject2.x01 > fitobject2.x02
        Fitted_D = [fitobject2.x01, fitobject2.x02];
        GaussFit1 = Gauss(fitobject2.sig1, fitobject2.x01, fitobject2.A1, Xfit);
        GaussFit2 = Gauss(fitobject2.sig2, fitobject2.x02, fitobject2.A2, Xfit);
    else
        Fitted_D = [fitobject2.x02, fitobject2.x01];
        GaussFit1 = Gauss(fitobject2.sig2, fitobject2.x02, fitobject2.A2, Xfit);
        GaussFit2 = Gauss(fitobject2.sig1, fitobject2.x01, fitobject2.A1, Xfit);
    end
    Yfit = Gauss2(fitobject2.sig1, fitobject2.x01, fitobject2.A1, fitobject2.sig2, fitobject2.x02, fitobject2.A2, Xfit);
    if UseLog==1; Fitted_D=10.^Fitted_D; end
    
    Idx_1 = UsedDapp(:,1) > x_Inter;
%     MSD_all_1 = MSD_all(Idx_1);
%     Traj_1 = Traj(Idx_1);
    Idx_2 = UsedDapp(:,1) < x_Inter;
%     MSD_all_2 = MSD_all(Idx_2);
%     Traj_2 = Traj(Idx_2);
%     MSD_1 = zeros(Lmax, 3); % MSD_1 corresponds to the distribution with the highest diffusion coeff
%     MSD_2 = zeros(Lmax, 3); % MSD_2 corresponds to the distribution with the lowest diffusion coeff
    
    UsedDapp_1 = UsedDapp(Idx_1,1);
    UsedDapp_2 = UsedDapp(Idx_2,1);
    F_mobile = round(100*size(UsedDapp_1,1)/(size(UsedDapp_1,1)+size(UsedDapp_2,1)));
    disp(strcat('The mobile fraction is equal to f=',num2str(F_mobile), '%'));
    
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
    plot(Xfit, GaussFit1, '-r', 'LineWidth',1)
    plot(Xfit, GaussFit2, '-b', 'LineWidth',1)
    plot(Xfit, Yfit, '--k', 'LineWidth',1)
    
    ax.FontSize = FontSize;
    axis square
    box on
    xlabel(x_label)
    ylabel('Fraction of tracks (%)')   
end

%3 Gaussians
if NbrFit == 3
    
    % Fit the three distributions with a Gaussian function and plot the
    % three fits as well as the sum of the three on the graph. 
    % ------------------------------------------------
    
    warnh = warndlg('Click on the approximate positions of the three maxima.','');
    uiwait(warnh)
    [X,Y] = ginput(3);
    
    X1 = X(1);
    A1 = 100*Y(1);
    X2 = X(2);
    A2 = 100*Y(2);
    X3 = X(3);
    A3 = 100*Y(3);
    
    warnh = warndlg('Click on the approximate positions of the FWHM for the first Gaussian.','');
    uiwait(warnh)
    [X,~] = ginput(2);
    sig1 = abs(X(2)-X(1))/2.35;
    
    warnh = warndlg('Click on the approximate positions of the FWHM for the second Gaussian.','');
    uiwait(warnh)
    [X,~] = ginput(2);
    sig2 = abs(X(2)-X(1))/2.35;
    
    warnh = warndlg('Click on the approximate positions of the FWHM for the third Gaussian.','');
    uiwait(warnh)
    [X,~] = ginput(2);
    sig3 = abs(X(2)-X(1))/2.35;
    
    Gauss3 = fittype( @(sig1,x01,A1,sig2,x02,A2,sig3,x03,A3,x) A1.*(exp(-((x - x01).^2)./(2*sig1^2))) ...
        + A2.*(exp(-((x - x02).^2)./(2*sig2^2))) + A3.*(exp(-((x - x03).^2)./(2*sig3^2))) );
    [fitobject3,gof3] = fit(double(BinCenters)', double(N)', Gauss3, 'start', [sig1,X1,A1,sig2,X2,A2,sig3,X3,A3]);
    disp(strcat('for the fit, R²=', num2str(100*gof3.rsquare), '%'))
    
    % Find the intersections points between the three fits (estimation)
    % -------------------------------------------------------------
    
    Xfit = min(BinCenters) : 0.01 : max(BinCenters);
    GaussFit1 = Gauss(fitobject3.sig1, fitobject3.x01, fitobject3.A1, Xfit);
    GaussFit2 = Gauss(fitobject3.sig2, fitobject3.x02, fitobject3.A2, Xfit);
    GaussFit3 = Gauss(fitobject3.sig3, fitobject3.x03, fitobject3.A3, Xfit);
    
    %between first and second
    X0=[fitobject3.x01, fitobject3.x02, fitobject3.x03];
    A=[fitobject3.A1, fitobject3.A2, fitobject3.A3];
    SIG=[fitobject3.sig1, fitobject3.sig2, fitobject3.sig3];
    
    [~,sorted_X0_index]=sort(X0); 
    X0=X0(sorted_X0_index);
    A=A(sorted_X0_index);
    SIG=SIG(sorted_X0_index);
    
    
    Idx_Bin = find(Xfit >= X0(1) & Xfit <= X0(2)); 
    Diff = abs(GaussFit1(Idx_Bin) - GaussFit2(Idx_Bin));
    [~,Idx] = min(Diff);
    x_Inter_1 = Xfit(Idx_Bin(Idx));
    
    Idx_Bin = find(Xfit >= X0(2) & Xfit <= X0(3)); 
    Diff = abs(GaussFit1(Idx_Bin) - GaussFit2(Idx_Bin));
    [~,Idx] = min(Diff);
    x_Inter_2 = Xfit(Idx_Bin(Idx));
    
    % Analyse separatly the three distributions and return the fraction of
    % events belonging to each of the three distributions. It also splits the
    % trajectories and MSD between the three populations.
    % -------------------------------------------------
       
    
    
    GaussFit1 = Gauss(SIG(1), X0(1), A(1), Xfit);
    GaussFit2 = Gauss(SIG(2), X0(2), A(2), Xfit);
    GaussFit3 = Gauss(SIG(3), X0(3), A(3), Xfit);
        
    Yfit = Gauss3(SIG(1), X0(1), A(1), SIG(2), X0(2), A(2), SIG(3), X0(3), A(3), Xfit);
    Fitted_D = X0;
    if UseLog==1; Fitted_D=10.^Fitted_D; end
    
    Idx_1 = UsedDapp(:,1) < x_Inter_1;
%     MSD_all_1 = MSD_all(Idx_1);
%     Traj_1 = Traj(Idx_1);
    
    Idx_2 = UsedDapp(:,1) >= x_Inter_1 & UsedDapp(:,1) < x_Inter_2;
%     MSD_all_2 = MSD_all(Idx_2);
%     Traj_2 = Traj(Idx_2);
    
    Idx_3 = UsedDapp(:,1) >= x_Inter_2;
%     MSD_all_3 = MSD_all(Idx_3);
%     Traj_3 = Traj(Idx_3);
    
%     MSD_1 = zeros(Lmax, 3); % MSD_1 corresponds to the distribution with the highest diffusion coeff
%     MSD_2 = zeros(Lmax, 3); % MSD_2 corresponds to the distribution with the lowest diffusion coeff
%     MSD_3 = zeros(Lmax, 3); % MSD_2 corresponds to the distribution with the lowest diffusion coeff
    
    UsedDapp_1 = UsedDapp(Idx_1,1);
    UsedDapp_2 = UsedDapp(Idx_2,1);
    UsedDapp_3 = UsedDapp(Idx_3,1);
    
    F_immobile = round(100*size(UsedDapp_1,1)/(size(UsedDapp_1,1)+size(UsedDapp_2,1)+size(UsedDapp_3,1)));
    disp(strcat('The immobile fraction is equal to f=',num2str(F_immobile), '%'));
    F_medium_mobile = round(100*size(UsedDapp_2,1)/(size(UsedDapp_1,1)+size(UsedDapp_2,1)+size(UsedDapp_3,1)));
    disp(strcat('The immobile fraction is equal to f=',num2str(F_medium_mobile), '%'));
    F_mobile = round(100*size(UsedDapp_3,1)/(size(UsedDapp_1,1)+size(UsedDapp_2,1)+size(UsedDapp_3,1)));
    disp(strcat('The immobile fraction is equal to f=',num2str(F_mobile), '%'));
    
    [N1,BinCenters] = hist(UsedDapp_1(:,1), BinCenters);
    [N2,BinCenters] = hist(UsedDapp_2(:,1), BinCenters);
    [N3,BinCenters] = hist(UsedDapp_3(:,1), BinCenters);
    
    N1_norm = 100*N1/sum(N1+N2+N3);
    N2_norm = 100*N2/sum(N1+N2+N3);
    N3_norm = 100*N3/sum(N1+N2+N3);
    
    % replot the two distributions, showing now the two separate
    % populations
    % -----------
    
    figure(hPlot)
    ax = gca;
    hold off
    cla
    
    b = bar(BinCenters', cat(2, N1_norm', N2_norm', N3_norm'), 'stacked');
    b(1).FaceColor = [1 0.5 0];
    b(2).FaceColor = [0 0.4 1];
    b(3).FaceColor = [1 0 1];
    
    hold on
    plot(Xfit, GaussFit1, '-r', 'LineWidth',1)
    plot(Xfit, GaussFit2, '-b', 'LineWidth',1)
    plot(Xfit, GaussFit3, '-g', 'LineWidth',1)
    plot(Xfit, Yfit, '--k', 'LineWidth',1)
    
    ax.FontSize = FontSize;
    axis square
    box on
    xlabel(x_label)
    ylabel('Fraction of tracks (%)')   
end
    
     if NbrFit==1          
         NewLine = sprintf('D = %.2f µm²/s', Fitted_D(1,:));
         Title = NewLine;
     elseif NbrFit==2          
         NewLine = sprintf('D_1 = %.2f µm²/s ; D_2 = %.2f µm²/s', Fitted_D(1,:));
         NewLine2 = sprintf('Mobile fraction : %d%%', F_mobile); 
         Title = char(NewLine, NewLine2);
     elseif NbrFit==3          
         NewLine = sprintf('D_1 = %.2f µm²/s ; D_2 = %.2f µm²/s ; D_3 = %.2f µm²/s', Fitted_D(1,:));
         NewLine2 = sprintf('Immobile fraction : %d%% ; Most mobile fraction : %d%%', F_immobile, F_mobile); 
         Title = char(NewLine, NewLine2);
     end
    
    title(Title);
    
% Save the plot
% -------------

% if DiffCalculationMethod == 1
%     export_fig(hPlot, 'Diffusion_distribution_Fit_Method.png');
%     export_fig(hPlot, 'Diffusion_distribution_Fit_Method.pdf', '-pdf');
% else
%     export_fig(hPlot, 'Diffusion_distribution_Weighted_Average_Method.png');
%     export_fig(hPlot, 'Diffusion_distribution_Weighted_Average_Method.pdf', '-pdf');
% end
    

