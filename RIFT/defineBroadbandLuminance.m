function S = defineBroadbandLuminance( t, f_c )
% Generates the luminance time series signal for broadband freq-tagging.
% See Zhigalov & Jensen (2020) Hum Brain Mapp. for the exact equation.
% 
% INPUT:
%   t  : vector, time variable.
%   f_c: scalar, frequency of the carrier signal.
% 
% OUTPUT:
%   S  : vector, luminance signal going from 0 to 1.
% 
% JY (Sep, 2022)

SIGMA = @(ii,t,r) sin(2*pi*t*(ii+r/4+1)); %sub-function

S = sin(2*pi*t*f_c + SIGMA(1,t,rand()) + SIGMA(2,t,rand()) + SIGMA(3,t,rand()) ); %this goes from -1 to 1
S = ( S - min(S) ) ./ 2; %rescales S so that it goes from 0 to 1.

end