% InitPert    EDIFICE: Initialise reference perturbation field
%
% [pert]  =  InitPert(FE,n,sym)
%
%   Initialise reference perturbation field to add noise to yield stress
%   for triggering localized deformation. Input langrangian particle
%   struct LP as set up by InitLP, finite element mesh struct FE as set
%   up by InitFE. Input n is the number of smoothing steps applied to white
%   noise; sym = 'ON' initialises a perturbation field that is symmetric
%   around a vertical line down the middle of the domain.
%
%   created   20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function  [pert]  =  InitPert(FE,n,sym)

rng(5,'twister');

%*****  produce random perturbations  *************************************

a  =  randn(FE.NIP,1)*0.3;
a  =  SmoothField(a,1,n,FE,'IP');

if strcmp(sym,'ON')
    a(FE.MapIP(:,1:(FE.nxIP-1)/2)) = fliplr(a(FE.MapIP(:,(FE.nxIP-1)/2+2:end)));
end

pert  =  a./max(abs(a));

