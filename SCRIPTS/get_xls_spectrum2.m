function spectrum = get_xls_spectrum2(xls_file, xls_sheet, xls_range, lambda)
% NAME:
%	GET_XLS_SPECTRUM
%
% PURPOSE:
%	Read a spectrum from Excel using readtable
%
% CATEGORY:
%	Signal, image processing.
%
% CALLING SEQUENCE:
%	spectrum = get_xls_spectrum(xls_file, xls_str, lambda);
%
% INPUTS:
%	xls_file: name of the excel file
%   xls_sheet: number containing the excel cell sheet to read
%   xls_range: string containing the excel cell range to read
%   lambda: normalize at this wavelength, if present. If not, normalize at
%           maximum
% OUTPUTS:
%	spectrum: the normalized spectrum
%
% COMMON BLOCKS:
%	None.
%
% SIDE EFFECTS:
%	None.
%
% RESTRICTIONS:
%	None.
%
% MODIFICATION HISTORY:
%	D.Bourgeois, March 2021.

if ~isempty(xls_range)
    spectrum = readtable(xls_file,'Sheet', xls_sheet,'Range', xls_range, 'PreserveVariableNames',true);
    spectrum = table2array(spectrum);
    if nargin == 4  %Normalize the spectrum at lambda
        tmp_var=abs(spectrum(:,1)-lambda);
        w_lambda=find(tmp_var == min(tmp_var));
        w_lambda=w_lambda(1); %index of the position where to normalize
        if spectrum(w_lambda,2)<=0
            disp(['Spectrum too weak at: ',num2str(lambda), 'nm']);
            disp('Cannot normalize: ');
            return;
        end
        spectrum(:,2)=spectrum(:,2)/spectrum(w_lambda,2);
        disp(['Spectrum normalized to 1 at: ',num2str(lambda), 'nm']);
    else %Normalize the spectrum at maximum
        tmp_var=abs(spectrum(:,2)-max(spectrum(:,2)));
        w_lambda=find(tmp_var == min(tmp_var));
        w_lambda=w_lambda(1); %index of the position where to normalize
        spectrum(:,2)=spectrum(:,2)/spectrum(w_lambda,2);
        disp(['Spectrum normalized to 1 at: ',num2str(spectrum(w_lambda,1)), 'nm']);
    end
    % Set to 0 negative values
    spectrum(spectrum(:,2)<0,2)=0;
    
    % Remove eventual nan values
    spectrum(isnan(spectrum(:,1)),:)=[];
    spectrum(isnan(spectrum(:,2)),:)=[];
else
    spectrum=0;
end
end

