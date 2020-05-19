% InitEN    FEMSTOKES: Set energy initial condition
%
% [CTX] = InitEN(CTX)
%
%   Function initializes the energy conservation solution variable for the
%   Stokes flow problem consisting of temperature T
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller


function  [CTX]  =  InitEN(CTX)


%***  prepare structures and variables
SL    =  CTX.SL;
FE    =  CTX.FE;
INIT  =  CTX.INIT;
BC    =  CTX.BC;
LP    =  CTX.LP;

%***  initialise with constant T = T0 if no energy coupling
if strcmp(SL.EN,'OFF'); INIT.TMode = 'constant'; end;

%***  prepare initial perturbation
pert  =  INIT.CosPert.*cos(INIT.NPert*2*pi.*FE.CoordQ1(:,1)./FE.W)  ...
      +  INIT.SinPert.*sin(INIT.NPert*2*pi.*FE.CoordQ1(:,1)./FE.W);

        
%*****  initialize temperature on particles  ******************************

SL.T  =  zeros(FE.NQ1,1);

%***  T = T0 with boundary layers towards TopTval and BotTval
if     strcmp(INIT.TMode(1:6),'blay') % 'blayer'
    
    dTbot  =  INIT.T0 - BC.BotTval;
    dTtop  =  INIT.T0 - BC.TopTval;
    
    SL.T   =  INIT.T0 - dTbot.*exp(      -FE.CoordQ1(:,2) ./INIT.Tthick) ...
                      - dTtop.*exp(-(FE.D-FE.CoordQ1(:,2))./INIT.Tthick);
    
    
%***  constant initial T = T0
elseif strcmp(INIT.TMode(1:4),'cons') % 'constant'
    
    SL.T(:)  =  INIT.T0;
    
    
%***  T-gradient from TopTval to BotTval
elseif strcmp(INIT.TMode(1:4),'grad') % 'gradient'
    
    SL.T  =  BC.TopTval - (FE.CoordQ1(:,2)./FE.D+pert) * (BC.TopTval-BC.BotTval);
    
    
%***  lithosphere-like piece-wise gradient with kink at LAB
elseif strcmp(INIT.TMode(1:4),'lith') % 'lithosphere'
    
    Tm    =  BC.BotTval + (FE.D-FE.CoordQ1(:,2)) * (INIT.T0-BC.BotTval)./(FE.D-INIT.LAB-pert);
    Tl    =  BC.TopTval +       FE.CoordQ1(:,2)  * (INIT.T0-BC.TopTval)./(     INIT.LAB+pert);
    SL.T  =  min(Tl,Tm);
    
    
%***  halfspace cooling model with constant age = TAge
elseif strcmp(INIT.TMode(1:4),'half') % 'halfspace'
    
    TAge   =  INIT.TAge;
    kappa  =  mean(CTX.PHYS.kS./CTX.MP.Rho(1)./CTX.PHYS.CpS);
    SL.T   =  BC.TopTval + (INIT.T0-BC.TopTval) .* erf(FE.CoordQ1(:,2)./(2.*sqrt(kappa.*TAge)));
    
%***  hot conduit through domain
elseif strcmp(INIT.TMode(1:4),'cond') % 'conduit'
    
    T0    =  BC.TopTval - (FE.CoordQ1(:,2)./FE.D+pert) * (BC.TopTval-BC.BotTval);
    cndt  =  (1 - tanh((FE.CoordQ1(:,1)-INIT.TXLoc-INIT.TWidth/2)/(INIT.TWidth/5)))/2.*(1+tanh((FE.CoordQ1(:,1)-INIT.TXLoc+INIT.TWidth/2)/(INIT.TWidth/5)))/2;
    SL.T  =  T0 + cndt.*(INIT.T0-T0);

%***  idealised lava lake bed
elseif strcmp(INIT.TMode(1:4),'lava') % 'lavalake'
    
    T0    =  BC.TopTval - (FE.CoordQ1(:,2)./FE.D+pert) * (BC.TopTval-BC.BotTval);
    width =  max(INIT.TWidth, FE.W-0.5 - (FE.W-0.5-INIT.TWidth)./INIT.THeight.*FE.CoordQ1(:,2));
    lake  =  FE.CoordQ1(:,1) >= INIT.TXLoc-width/2 & FE.CoordQ1(:,1) <= INIT.TXLoc+width/2;
    SL.T  =  T0 + lake.*(INIT.T0-T0);
    
end


%***  add a gaussian T-pulse on top of the T-initial condition
if strcmp(INIT.TMode(end-4:end),'pulse')
    Tpulse  =  INIT.TAmpl .* exp(-(FE.CoordQ1(:,1)-INIT.TXLoc).^2 ./ INIT.TWidth^2) ...
                          .* exp(-(FE.CoordQ1(:,2)-INIT.TZLoc).^2 ./ INIT.THeight^2);
else
    Tpulse  =  0;
end

Trand  =  INIT.PertTemp.*PLPQ1(LP.pert,LP,FE,LP.LP2FEmethod);

SL.T   =  SL.T + Tpulse + Trand;

SL.T  =  SmoothField(SL.T,0.125,round(CTX.FE.nxEl^2/500),FE,'Q1');


%*****  interpolate temperature to particles  *****************************

if strcmp(CTX.SL.Advection,'LGR')
    LP.T  =  PQ1LP(SL.T,LP,FE,LP.FE2LPmethod);
end

%*****  store boundary T for constant T BC  *******************************

if strcmp(BC.TopIsotherm,'CC')    % top    boundary retains initial distribution
    BC.TopTval   = SL.T(FE.MapQ1(1  ,:)).';
end
if strcmp(BC.BotIsotherm,'CC')    % bottom boundary retains initial distribution
    BC.BotTval   = SL.T(FE.MapQ1(end,:)).';
end
if strcmp(BC.LeftIsotherm,'CC')   % left   boundary retains initial distribution
    BC.LeftTval  = SL.T(FE.MapQ1(:,1  ));
end
if strcmp(BC.RightIsotherm,'CC')  % right  boundary retains initial distribution
    BC.RightTval = SL.T(FE.MapQ1(:,end));
end


%***  return structures
CTX.SL  =  SL;
CTX.BC  =  BC;
CTX.LP  =  LP;

end



