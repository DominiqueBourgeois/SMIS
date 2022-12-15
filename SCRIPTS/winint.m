function x = winint(x,int,m,v);

% WININT  Returns a normalized Window of Values in Intervall
%
% Y = WININT( X , Intervall , Mode , Deviation )
%
% Mode    = Mode of window: {'triangle'} | 'cosine' |  'gauss'  |  'exp'   
%
%           gauss: exp( -(pi*x)^2 / 2 );  exp: exp(-x);
%
% Mode    = positive Numeric for Power of Window: ( -1 .. 1 )
%             Mode == 0    ==>  Impuls
%             Mode == 1    ==>  Triangle
%             Mode == Inf  ==>  Linear (constant)
% 
% Deviation = shift of Window:
%              -1 left    ¯¯\_          1 ---> 0
%               0 center  _/\_   0 ---> 1 ---> 0  (default)
%               1 right   _/¯¯   0 ---> 1
%
%
%  DEV   Y(X<=INT(1))  Y(X>=INT(2))  Y(X==mean(INT))
% ---------------------------------------------------
%   -1   1             0
%    0   0             0             1
%    1   0             1
%
% In case of Mode "exp" the Values nears asymtotic to ZERO!
%
%----------------------------------------------------
% Example:
%
%  x = ( 0 : 100 );  int = [ 30  70 ];
%  c = {'r--' 'k-' 'g--'};
%  for mode = { 'triangle' 'cosine' 'gauss'  2  0.5  0 inf }
%      figure, hold on, box on
%      xlim([0 100]), ylim([0 1]+0.05*[-1 1])
%      for v = [ 0  -1   1 ]
%          plot(x,winint(x,int,mode{1},v),c{v+2});
%      end
%      if ischar(mode{1})
%         title(mode{1})
%      else
%         title(sprintf('Potenz: %.3g',mode{1}));
%      end
%  end
%
%               
%----------------------------------------------------
%
% see also: WINDOW, PLATEAU
%
       
Nin = nargin;

if Nin < 2
   error('Not enough InputArguments.');
end

msg = cell(0,1);

%---------------------------------------------------
% Check X

if ~isnumeric(x)
    msg = cat(1,msg,{'X must be numeric.'});
elseif ~strcmp(class(x),'double')
    x = double(x);
end

%---------------------------------------------------
% Check Intervall

if ~( isnumeric(int) & ( prod(size(int)) == 2 ) )
    msg = cat(1,msg,{'Intervall must be a 2 Element numeric.'});
elseif ~strcmp(class(int),'double')
    int = double(int);
end


%---------------------------------------------------
% Check Mode

e = 1;

if Nin < 3

  m = 'l';

else

  ok = ( ischar(m) & ( prod(size(m)) == size(m,2) ) & ~isempty(m) );

  if ok
     m = lower(m(1));
  else
     ok = ( isnumeric(m)  &  ( prod(size(m)) == 1 ) );
     if ok
        ok = ( ~isnan(m) & ( m >= 0 ) );
     end
     if ok
        e = m;
        m = 'p';
     end
  end

  if ~ok
      msg = cat( 1 , msg , {'Mode must be a String or positive Numeric.'} );
  end

end

%---------------------------------------------------
% Check Deviation

if Nin < 4

  v = 0;

else

  if ~( isnumeric(v)  &  ( prod(size(v)) == 1 ) );
     msg = cat( 1 , msg , {'Potenz must be a finite Numeric.'} );
  else
     v = sign(double(v));
  end 

end

%---------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%*****************************************************************

int = int( [ 1  2 ] + [ 1 -1 ] * ( int(2) < int(1) ) );

%-----------------------------------------------------
% Offset X0
%
% v  = [  -1        0        1   ];
% x0 = [ int(1) mean(int) int(2) ]


x0 = [ 1  2 ];
if ~( v == 0 )
     ii  = ( 3 - v ) / 2;
  x0(ii) = x0(ii) + v;
end

x0 = int(x0);

x0 = ( x0(1) + x0(2) ) / 2;

%-----------------------------------------------------
% Check for ZERO-Intervall

nn = find(isnan(x));

if int(1) == int(2)

   if v == 0
      x = ( x == x0 );
   else
      x = ( v == -1 ) + v * ( x > x0 );
   end

   if ~isempty(nn)
      x(nn) = NaN;
   end

   return   

end

x = max(x,int(1));

if 1 %%% ~any( m == 'ge' )
   x = min( x , int(2) );
end

%-----------------------------------------------------
% Scale with Period
%
% v = [     -1           0           1     ]
% x = { [  0   1 ]  [ -1   1 ]  [ -1   0 ] }

x = ( 1 + (v==0) ) * ( x - x0 ) / ( int(2) - int(1) );


%-----------------------------------------------------

switch m

  %-------------------------------------------------
  % Cosine

  case 'c'

    x = ( 1 + cos( pi * x ) ) / 2;

  %-------------------------------------------------
  % Gauss

  case 'g'

    x = exp( (-1) * (pi*x).^2 / 2 );

    b = exp( (-1) * (pi*1).^2 / 2 );

    x = ( x - b ) / ( 1 - b );

%%%    x = ( x - 1 ) / ( 1 - exp( (-1) * (pi).^2 / 2 ) ) + 1;

  %-------------------------------------------------
  % Exponent

  case 'e'

    x = exp( (-1) * abs(x) );

  %-------------------------------------------------
  % Triangle | Potenz

  otherwise

    if isinf(e)

       x = ( abs(x) < 1 );

    else
 
       if e == 0
          x = double( x == 0 );
       elseif e == 1
          x = 1 - abs(x);
       else    
          x = 1 - abs(x) .^ e;
       end

    end

end

if ~isempty(nn)
    x(nn) = NaN;
end
