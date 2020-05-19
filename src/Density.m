
% Density    EDIFICE: Update material density
%
% [CTX] = Density(MP,CTX)
%
%   Function updates T-dependent density according to latest solution guess. 
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function   [MP]  =  Density(MP,CTX)


%*****  get densities from material types  ********************************

MP.Rho    =  CTX.PROP.Rho(CTX.MP.Mat);


end





