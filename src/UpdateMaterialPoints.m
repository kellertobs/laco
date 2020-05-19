
% UpdateMaterialProps    EDIFICE: Update material properties according to current solution
%
% [CTX] = UpdateMaterialProps(CTX)
%
%   Function updates material properties and auxiliary fields according to
%   latest solution guess and parameters and options provided in the CTX
%   struct.
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function [CTX]  =  UpdateMaterialPoints(CTX)

MP  =  CTX.MP;


%*****  Update density  ***************************************************

MP.Rho    =  CTX.PROP.Rho(CTX.MP.Mat);


%*****  Update stress / strain rates  *************************************

MP  =  StressStrainr(MP,CTX);


%*****  Update rheology  **************************************************

MP  =  Rheology(MP,CTX);


CTX.MP  =  MP;

end







