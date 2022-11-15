function  sm_par = get_sm_filter_profiles(im_par,sm_par)

% NAME:
%	get_filter_profiles
%
% PURPOSE:
%	Get filter profiles for specific fluorescence emission spectra so that
%	they can be multiplied point by point
%
% INPUTS:
%   im_par: the SMIS imaging parameters
%	sm_par: the SM parameters
%
% OUTPUTS:
%	sm_par: the updated SM parameters for .filter_profiles
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022


plot_profiles=0; % To see the profiles

em_f1 = im_par.filters.ch1_filter; % Get filter 1

if im_par.two_channel==1
    em_f2 = im_par.filters.ch2_filter; % Get filter 2
end

filter_profiles(1:sm_par.n_fluorescent_states)=struct;

for i=1:sm_par.n_fluorescent_states
    em_spectrum=sm_par.spectral_data.em_spectra(i).s; % Reference emission spectrum
    l=em_spectrum(:,1); % The wavelengths

    % Emission filter 1
    if im_par.filters.use_exp_ch1_filter==1
        x_ch=im_par.filters.exp_ch1_filter(1,:);
        y_ch=im_par.filters.exp_ch1_filter(2,:);
    else
        [x_ch,y_ch,~] = get_filter_profile(em_f1);
    end

    if min(l)<min(x_ch) && max(l)>max(x_ch) % General case
        xf=[min(l), x_ch, max(l)];
        yf=[0, y_ch, 0];
    elseif min(l)>=min(x_ch) && max(l)>max(x_ch)
        xf=[x_ch, max(l)];
        yf=[y_ch, 0];
    elseif min(l)<min(x_ch) && max(l)<=max(x_ch)
        xf=[min(l), x_ch];
        yf=[0, y_ch];
    else
        xf=x_ch;
        yf=y_ch;
    end
    filter_profiles(i).ch1=interp1(xf,yf,l, "linear"); % Get the filter values at the lambdas of the spectrum

    % Emission filter 2
    if im_par.two_channel==1

        if im_par.filters.use_exp_ch2_filter==1
            x_ch=im_par.filters.exp_ch2_filter(1,:);
            y_ch=im_par.filters.exp_ch2_filter(2,:);
        else
            [x_ch,y_ch,~] = get_filter_profile(em_f2);
        end

        if min(l)<min(x_ch) && max(l)>max(x_ch) % General case
            xf=[min(l), x_ch, max(l)];
            yf=[0, y_ch, 0];
        elseif min(l)>=min(x_ch) && max(l)>max(x_ch)
            xf=[x_ch, max(l)];
            yf=[y_ch, 0];
        elseif min(l)<min(x_ch) && max(l)<=max(x_ch)
            xf=[min(l), x_ch];
            yf=[0, y_ch];
        else
            xf=x_ch;
            yf=y_ch;
        end
        filter_profiles(i).ch2=interp1(xf,yf,l, "linear"); % Get the filter values at the lambdas of the spectrum

        if im_par.filters.add_dichroic==1

            if im_par.filters.use_exp_dichroic==1
                x_ch=im_par.filters.exp_dichroic_filter(1,:);
                y_ch=im_par.filters.exp_dichroic_filter(2,:);
                if min(l)<min(x_ch) && max(l)>max(x_ch) % General case
                    xf=[min(l), x_ch, max(l)];
                    yf=[0, y_ch, 0];
                elseif min(l)>=min(x_ch) && max(l)>max(x_ch)
                    xf=[x_ch, max(l)];
                    yf=[y_ch, 0];
                elseif min(l)<min(x_ch) && max(l)<=max(x_ch)
                    xf=[min(l), x_ch];
                    yf=[0, y_ch];
                else
                    xf=x_ch;
                    yf=y_ch;
                end
            else
                df=im_par.filters.dichroic_filter;
                if min(l)<(df(1)-df(2)/2) && max(l)>(df(1)+df(2)/2) % General case
                    xf=[min(l),df(1)-df(2)/2,df(1)+df(2)/2, max(l)];
                    yf=[df(3), df(3), df(4),df(4)];
                elseif min(l)<(df(1)-df(2)/2) && max(l)<=(df(1)+df(2)/2)
                    xf=[min(l),df(1)-df(2)/2,df(1)+df(2)/2];
                    yf=[df(3), df(3), df(4)];
                elseif min(l)>=(df(1)-df(2)/2) && max(l)>(df(1)+df(2)/2)
                    xf=[df(1)-df(2)/2,df(1)+df(2)/2,max(l)];
                    yf=[df(3), df(4), df(4)];
                else
                    xf=[df(1)-df(2)/2,df(1)+df(2)/2];
                    yf=[df(3), df(4)];
                end
            end
        end
        
        filter_profiles(i).df=interp1(xf,yf,l, "linear"); % Get the filter values at the lambdas of the spectrum
    end

    if plot_profiles==1
        s=em_spectrum(:,2); % the spectrum
        s=s/max(s);
        plot(l,s,'blue','LineWidth',2);
        hold on
        plot(l,filter_profiles(i).ch1,'green','LineWidth',2);
        if im_par.two_channel==1
            plot(l,filter_profiles(i).ch2,'red','LineWidth',2);
            if im_par.filters.add_dichroic==1
                plot(l,filter_profiles(i).df,'magenta','LineWidth',2);
            end
        end
        xlabel('Wavelength [nm]')
        ylabel('Normalized intesity [AU]')
        ylim([0,1.1])
        set(gca,'LineWidth',2,'FontSize',16)

    end
end

sm_par.filter_profiles=filter_profiles;

