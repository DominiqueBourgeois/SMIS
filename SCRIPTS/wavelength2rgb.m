function [c,w] = wavelength2rgb(n,scl)

% SPECTRA  RGB-Values of Visible ColorSpectra
%
% [ RGB , W ]  = SPECTRA( N )
%
% Returns the RGB-Values and WaveLengths for
%  an N-point ColorSpectra between 380nm <= W <= 780nm.
%
% [ RGB , W ] = SPECTRA( [ W0 W1 ] ); with  N = W1-W0+1
% [ RGB , W ] = SPECTRA( [ W0 W1 N ] )
%
% Returns the RGB-Values and WaveLengths for
%  an N-point ColorSpectra between W0 <= W <= W1.
%
% [ RGB , W ] = SPECTRA( W )
%
%  Returns the RGB-Values for the WaveLengths defined by W,
%    use SPECTRA(-W) for a single WaveLength,
%    use a Column-Vector for two or three elements in W.
%
% Near the Vision Limits (below 420nm and above 700nm)
% The Intensity falls by linear by 0.7, the complete
% is GAMMA-Adjusted by the power of 0.8.
% A second Input redefines that both Values:
%
%   SPECTRA( ... , [ IntensFall  Gamma ] )
%
%--------------------------------------------------------
%
% RGB VALUES FOR VISIBLE WAVELENGTHS
%  by Dan Bruton (astro@tamu.edu)
%
%  http://www.physics.sfasu.edu/astro/color.html
%
%--------------------------------------------------------
%
% see also: BBODY
%

if nargin < 1
   n = [];
end

lim = [ 380  780 ];

lin = 0.7;  % LET THE INTENSITY FALL OFF NEAR THE VISION LIMITS
pot = 0.8;  % GAMMA ADJUST

w = [];

p = prod(size(n));

if     p == 2
   lim = n;
     n = [];
elseif p == 3
   lim = n([1 2]);
     n = n(3);
elseif p > 3
     w = n(:);
     n = [];
end

if ~isempty(n)
    ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
    if ok
       if n < 0
          w = -n;
          n = 1;
       else
          ok = ( ( n > 0 ) & ( mod(n,1) == 0 ) );
       end
    end
    if ~ok
        error('N must be an Integer larger Zero.');
    end
end


if nargin == 2
   if ~( isnumeric(scl) & ( prod(size(scl)) == 2 ) )
       error('Second Input must be a 2-element Numeric.');
   end
   lin = scl(1);
   pot = scl(2);
end

if isempty(w)
   if isempty(n)
     n = ceil(abs(diff(lim)));
     n = max(1,n);
   end
   w = linspace(lim(1),lim(2),n)';
end

n = size(w,1);

c = zeros(n,3);

c(:,1) = winint(w,[380 440],'tri',-1) + ...
         winint(w,[510 580],'tri', 1);

c(:,2) = 1 - ( winint(w,[440 490],'tri',-1) + ...
               winint(w,[580 645],'tri', 1) ) ;

c(:,3) = winint(w,[490 510],'tri',-1);

fak = 1 - lin * ( winint(w,[380 420],'tri',-1) + ...
                  winint(w,[700 780],'tri', 1) );

c = ( c .* fak(:,[1 1 1]) ) .^ pot;



