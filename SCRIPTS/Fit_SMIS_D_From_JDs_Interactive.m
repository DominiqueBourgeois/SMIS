function [NbrFit, Fitted_D, Xfit, Yfit, Fitted_F] = Fit_SMIS_D_From_JDs_Interactive(N, BinEdges, UsedDapp, tracklength, FontSize, x_label, hPlot)

% PURPOSE:
%	Routine to interactively fit diffusion coefficients histograms calculated from jump distances
%
% INPUTS:
%   N: Normalized histogram values
%   BinEdges: Histogram bins edges
%   UsedDapp: The used apparent diffusion coefficients
%   Tracklength: The track length used to calculate the histogram
%   FontSize: The font size
%   x_label: The X label
%   hPlot: Handle to the used figure
%
% OUTPUTS:
%   NbrFit: Number of diffusion coefficients fitted
%   Fitted_D: Values of fitted diffusion coefficients
%   Xfit: The fitted X data
%   Yfit: The fitted Y data
%   Fitted_F: Fraction of mobile molecules (1 if 1 diffusion regimes, [immobile,mobile] fraction if 2, [immobile,mobile,most mobile] fractions if 3)
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022. Compatibility with SMIS 1.3. 


NbrFit = questdlg('How many diffusion regimes ?', 'Fit', '1', '2','3','1');
NbrFit = str2double(NbrFit);

BinWidth=BinEdges-circshift(BinEdges,[0,1]);
BinCenters=BinEdges(1:end-1)+BinWidth(2:end)/2;

N=100*N; % Go in percentage

n=tracklength;

% Define function
F = fittype( @(D,A,x) A.*1/factorial(n-1).*(n./D).^n.*x.^(n-1).*exp(-n.*x/D));

%Only one diffusion coefficient
if NbrFit == 1
    
    [A0,Idx] = max(N); % Estimate the starting parameters for the first fit
    A0 = double(A0);
    D0 = double(BinCenters(Idx));
    X2Fit=double(BinCenters)';
    Y2Fit=double(N)';
    
    LowerBounds=[0 0];
    UpperBounds=[inf inf];
   
    [fitobject,gof] = fit(X2Fit, Y2Fit, F, 'start', [D0,A0], 'Lower', LowerBounds,'Upper',UpperBounds);
    
    disp(strcat('For the fit, R²=', num2str(100*gof.rsquare), '%'))
    
    Xfit = min(BinCenters) : 0.01 : max(BinCenters);
    Yfit = F(fitobject.D, fitobject.A, Xfit);
    Fitted_D = fitobject.D;
    
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

% 2 diffusion coefficients
if NbrFit == 2
    
    % Define function
    F2 = fittype( @(D1,A1,D2,A2,x) A1.*1/factorial(n-1).*(n./D1).^n.*x.^(n-1).*exp(-n.*x/D1) + ...
        A2.*1/factorial(n-1).*(n./D2).^n.*x.^(n-1).*exp(-n.*x/D2));
    
    warnh = warndlg('Click on the approximate positions of the two maxima.','');
    uiwait(warnh)
    [X,Y] = ginput(2);
    
    D1 = X(1);
    A1 = Y(1);
    D2 = X(2);
    A2 = Y(2);
    
    X2Fit=double(BinCenters)';
    Y2Fit=double(N)';
    
    LowerBounds=[0 0 0 0];
    UpperBounds=[inf inf inf inf];

    [fitobject2,gof2] = fit(X2Fit, Y2Fit, F2, 'start', [D1,A1,D2,A2],'Lower', LowerBounds,'Upper',UpperBounds);
    disp(strcat('for the fit, R²=', num2str(100*gof2.rsquare), '%'))
    
    Xfit = min(BinCenters) : 0.01 : max(BinCenters);
    
    Yfit1 = F(fitobject2.D1, fitobject2.A1, Xfit);
    Yfit2 = F(fitobject2.D2, fitobject2.A2, Xfit);
    
    if fitobject2.D1 < fitobject2.D2
        Idx_Bin = find(Xfit >= fitobject2.D1 & Xfit <= fitobject2.D2);
    else
        Idx_Bin = find(Xfit >= fitobject2.D2 & Xfit <= fitobject2.D1);
    end
    
    Diff = abs(Yfit1(Idx_Bin) - Yfit2(Idx_Bin));
    [~,Idx] = min(Diff);
    x_Inter = Xfit(Idx_Bin(Idx));
    
    % Analyse separatly the two distributions and return the fraction of
    % events belonging to each of the two distributions. It also splits the
    % trajectories and MSD between the two populations.
    % -------------------------------------------------
    
    if fitobject2.D1 > fitobject2.D2
        Fitted_D = [fitobject2.D1, fitobject2.D2];
        Yfit1 = F(fitobject2.D1, fitobject2.A1, Xfit);
        Yfit2 = F(fitobject2.D2, fitobject2.A2, Xfit);
    else
        Fitted_D = [fitobject2.D2, fitobject2.D1];
        Yfit1 = F(fitobject2.D2, fitobject2.A2, Xfit);
        Yfit2 = F(fitobject2.D1, fitobject2.A1, Xfit);
    end
    
    Yfit = F2(fitobject2.D1, fitobject2.A1, fitobject2.D2, fitobject2.A2, Xfit);
       
    Idx_1 = UsedDapp(:,1) > x_Inter;
    
    Idx_2 = UsedDapp(:,1) < x_Inter;
    
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
    plot(Xfit, Yfit1, '-r', 'LineWidth',1)
    plot(Xfit, Yfit2, '-b', 'LineWidth',1)
    plot(Xfit, Yfit, '--k', 'LineWidth',1)
    
    ax.FontSize = FontSize;
    axis square
    box on
    xlabel(x_label)
    ylabel('Fraction of tracks (%)')
    
end

% 3 diffusion coefficients
if NbrFit == 3
    
    % Define function
    F3 = fittype( @(D1,A1,D2,A2,D3,A3,x) A1.*1/factorial(n-1).*(n./D1).^n.*x.^(n-1).*exp(-n.*x/D1) + ...
        A2.*1/factorial(n-1).*(n./D2).^n.*x.^(n-1).*exp(-n.*x/D2) + ...
        A3.*1/factorial(n-1).*(n./D3).^n.*x.^(n-1).*exp(-n.*x/D3));
    
    warnh = warndlg('Click on the approximate positions of the three maxima.','');
    uiwait(warnh)
    [X,Y] = ginput(3);
    
    D1 = X(1);
    A1 = 100*Y(1);
    D2 = X(2);
    A2 = 100*Y(2);
    D3 = X(3);
    A3 = 100*Y(3);
    
    X2Fit=double(BinCenters)';
    Y2Fit=double(N)';
    
    LowerBounds=[0 0 0 0 0 0];
    UpperBounds=[inf inf inf inf inf inf];
   
    [fitobject3,gof3] = fit(X2Fit, Y2Fit, F3, 'start', [D1,A1,D2,A2,D3,A3],'Lower', LowerBounds,'Upper',UpperBounds);
    disp(strcat('for the fit, R²=', num2str(100*gof3.rsquare), '%'))
    
    Xfit = min(BinCenters) : 0.01 : max(BinCenters);
    
    Yfit1 = F(fitobject3.D1, fitobject3.A1, Xfit);
    Yfit2 = F(fitobject3.D2, fitobject3.A2, Xfit);
%     Yfit3 = F(fitobject3.D3, fitobject3.A3, Xfit);
    
    
    D0=[fitobject3.D1, fitobject3.D2, fitobject3.D3];
    A=[fitobject3.A1, fitobject3.A2, fitobject3.A3];
    
    [~,sorted_D0_index]=sort(D0);
    D0=D0(sorted_D0_index);
    A=A(sorted_D0_index);
    
    Idx_Bin = find(Xfit >= D0(1) & Xfit <= D0(2));
    Diff = abs(Yfit1(Idx_Bin) - Yfit2(Idx_Bin));
    [~,Idx] = min(Diff);
    x_Inter_1 = Xfit(Idx_Bin(Idx));
    
    Idx_Bin = find(Xfit >= D0(2) & Xfit <= D0(3));
    Diff = abs(Yfit1(Idx_Bin) - Yfit2(Idx_Bin));
    [~,Idx] = min(Diff);
    x_Inter_2 = Xfit(Idx_Bin(Idx));
    
    
    % Analyse separatly the three distributions and return the fraction of
    % events belonging to each of the three distributions. It also splits the
    % trajectories and MSD between the three populations.
    % -------------------------------------------------
    
    Yfit1 = F(D0(1), A(1), Xfit);
    Yfit2 = F(D0(2), A(2), Xfit);
    Yfit3 = F(D0(3), A(3), Xfit);
    
    Yfit = F3(D0(1), A(1), D0(2), A(2), D0(3), A(3), Xfit);
    
    Fitted_D = D0;
       
    Idx_1 = UsedDapp(:,1) < x_Inter_1;
    
    Idx_2 = UsedDapp(:,1) >= x_Inter_1 & UsedDapp(:,1) < x_Inter_2;
    
    Idx_3 = UsedDapp(:,1) >= x_Inter_2;
    
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
    plot(Xfit, Yfit1, '-r', 'LineWidth',1)
    plot(Xfit, Yfit2, '-b', 'LineWidth',1)
    plot(Xfit, Yfit3, '-g', 'LineWidth',1)
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
    Fitted_F=1;
elseif NbrFit==2
    NewLine = sprintf('D_1 = %.2f µm²/s ; D_2 = %.2f µm²/s', Fitted_D(1,:));
    NewLine2 = sprintf('Mobile fraction : %d%%', F_mobile);
    Title = char(NewLine, NewLine2);
    Fitted_F=[1-0.01*F_mobile, 0.01*F_mobile];
elseif NbrFit==3
    NewLine = sprintf('D_1 = %.2f µm²/s ; D_2 = %.2f µm²/s ; D_3 = %.2f µm²/s', Fitted_D(1,:));
    NewLine2 = sprintf('Immobile fraction : %d%% ; Most mobile fraction : %d%%', F_immobile, F_mobile);
    Title = char(NewLine, NewLine2);
    Fitted_F=[0.01*F_immobile, 1-0.01*F_immobile-0.01*F_mobile, 0.01*F_mobile];
end


title(Title);

