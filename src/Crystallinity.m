

function  [CTX]  =  Crystallinity(CTX)

MP      =  CTX.MP;
MPo     =  CTX.MPo;
FE      =  CTX.FE;
M       =  CTX.MODEL;
dt      =  CTX.TIME.step;

kappa   =  M.k;
T0      =  mean(CTX.BC.TopTval);
tau     =  M.TopCoolTime;
delta   =  M.TopCoolDepth;
src     =  MP.Psi.*CTX.PHYS.DissFact./CTX.SL.RhoRef./CTX.PHYS.Cp ...
        -  (CTX.SL.T - T0) ./ max(2.*dt,tau) .* exp(-FE.CoordQ1(:,2)./delta)  ...
        +  MP.Gamma .* M.LatHeatdT;
srco    =  CTX.MPo.Psi.*CTX.PHYS.DissFact./CTX.SL.RhoRef./CTX.PHYS.Cp ...
        -  (CTX.SLo.T - T0) ./ max(2.*dt,tau) .* exp(-FE.CoordQ1(:,2)./delta)  ...
        +  CTX.MPo.Gamma .* M.LatHeatdT;
vel     =  [CTX.SL .U,CTX.SL .W];
velo    =  [CTX.SLo.U,CTX.SLo.W];
dv      =  zeros(size(vel));

[CTX.SL.T]  =  ExplicitAdvDiffSolver(CTX.SL.T,CTX.SLo.T,vel,velo,dv,dv,src,srco,kappa,dt,'Q1',CTX);

if M.IOFlowRate>0
    if FE.NU == FE.NQ2
        W   =  -PQ2Q1(CTX.SL.W,FE);
    else
        W   =  -CTX.SL.W;
    end
    ind  =  W(FE.MapQ1(end,:)) > 1e-3*M.IOFlowRate;
    CTX.SL.T(FE.MapQ1(end,ind))  =  CTX.MODEL.IOFlowTemp;
end

T         =  CTX.SL.theta_dt .* CTX.SL.T   + (1-CTX.SL.theta_dt) .* CTX.SLo.T;
chi       =  CTX.SL.theta_dt .* CTX.MP.chi + (1-CTX.SL.theta_dt) .* CTX.MPo.chi;
chieq     =  max(0,min(1,(T - M.Tliq)./(M.Tsol - M.Tliq))).^M.q;
Gamma     =  (chieq-chi)./dt;
MP.Gamma  =  0.5.*Gamma + 0.5.*MP.Gamma;

kappa   =  0;
src     =  MP.Gamma;
srco    =  MPo.Gamma;

DRho    =  M.RhoM - M.RhoX;
MP.DWC  =  - 2/9.*DRho.*CTX.PHYS.grav(2).*M.b0^2./MP.EtaVEP.*(1-MP.phi).^M.m;
vel     =  [CTX.SL .U,CTX.SL .W];
velo    =  [CTX.SLo.U,CTX.SLo.W];
dv      =  [0.*CTX.SL .U,MP .DWC];
dvo     =  [0.*CTX.SLo.U,MPo.DWC];

MP.chi  =  ExplicitAdvDiffSolver(MP.chi,CTX.MPo.chi,vel,velo,dv,dvo,src,srco,kappa,dt,'Q1',CTX);
MP.chi  =  max(0,min(1,MP.chi));

CTX.MP  =  MP;

end