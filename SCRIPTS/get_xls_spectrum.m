function spectrum = get_xls_spectrum(xls_file, xls_str, lambda)
% NAME:
%	GET_XLS_SPECTRUM
%
% PURPOSE:
%	Read a spectrum from Excel
%
% CATEGORY:
%	Signal, image processing.
%
% CALLING SEQUENCE:
%	spectrum = get_xls_spectrum(xls_file, xls_str, lambda);
%
% INPUTS:
%	xls_file: name of the excel file
%   xls_str: string containing the excel cell to read
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
%	D.Bourgeois, April 2011.

if ~isempty(xls_str)
    spectrum = xlsread(xls_file,1,xls_str);
    if nargin == 3  %Normalize the spectrum at lambda
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
else
    spectrum=0;
end
end

