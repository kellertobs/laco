

function  [L,RHS,BC]  =  FMApplyBC(L,RHS,X,BC,FE,CTX)

mapv   =  FE.MapV;
mapP   =  FE.MapP;

apply_topvx        =  0;
apply_topvz        =  0;
apply_botvx        =  0;
apply_botvz        =  0;
apply_leftvx       =  0;
apply_leftvz       =  0;
apply_rightvx      =  0;
apply_rightvz      =  0;
apply_Pfix         =  0;

if     strcmp(BC.TopBot,'ns')  % no slip
    apply_topvx  =  1;
    apply_topvz  =  1;
    apply_botvx  =  1;
    apply_botvz  =  1;
    
elseif strcmp(BC.TopBot,'fs')  % free slip
    apply_topvz  =  1;
    apply_botvz  =  1;
    
elseif strcmp(BC.TopBot,'bf')  % bottom flux
    apply_topvz  =  1;
    apply_botvz  =  1;
    
elseif strcmp(BC.TopBot,'ss')  % simple shear
    apply_topvx  =  1;
    apply_botvx  =  1;
    apply_topvz  =  1;
    apply_botvz  =  1;
end


if     strcmp(BC.Sides,'ns')  % no slip
    apply_leftvx  =  1;
    apply_leftvz  =  1;
    apply_rightvx =  1;
    apply_rightvz =  1;
    
elseif strcmp(BC.Sides,'fs')  % free slip
    apply_leftvx  =  1;
    apply_rightvx =  1;
    
elseif strcmp(BC.Sides,'rs')  % right side only free slip
    apply_leftvx  =  1;
    apply_leftvz  =  1;
    apply_rightvx =  1;
    
elseif strcmp(BC.Sides,'ls')  % left side only free slip
    apply_leftvx  =  1;
    apply_rightvx =  1;
    apply_rightvz =  1;
    
elseif strcmp(BC.Sides,'ss')  % simple shear
    apply_leftvx  =  1;
    apply_rightvx =  1;
    apply_leftvz  =  1;
    apply_rightvz =  1;
    
elseif strcmp(BC.Sides,'ff')  % free horizontal flux
    apply_leftvz  =  1;
    apply_rightvz =  1;  
end

if strcmp(BC.FreeSurface,'ON')
    apply_topvx    =  0;
    apply_topvz    =  0;
else
    apply_Pfix     =  1;
end

if strcmp(BC.Type,'ConstStr')
    if strcmp(BC.Sides,'ss')
        BC.UTopBot(1,:) =  - 0.5.*BC.BGStrainr.*FE.D;
        BC.UTopBot(2,:) =  + 0.5.*BC.BGStrainr.*FE.D;
    else
        BC.USides(1,:)  =  + 0.5.*BC.BGStrainr.*FE.W;
        BC.USides(2,:)  =  - 0.5.*BC.BGStrainr.*FE.W;
    end
end
    
if strcmp(BC.Sides,'ss')
    
    dv = BC.UTopBot(2) - BC.UTopBot(1);
        
    j  =  (1:FE.nzU)';
    U  =  BC.UTopBot(1) + repmat((j-1)*dv/FE.nzU,1,FE.nxU);
    W  =  zeros(size(SL.U));
    
else
    
    dv  =  BC.USides(2) - BC.USides(1);
    
    i   =  1:FE.nxU;
    U   =  BC.USides(1) + repmat((i-1)*dv/(FE.nxU-1),FE.nzU,1);
    
    j   =  (1:FE.nzU)';
    W   =  -repmat((j-1)*dv/FE.aspect/(FE.nzU-1),1,FE.nxU);
    
end
 
BC.UTopBot       =  zeros(2,FE.nxU);
BC.WTopBot       =  zeros(2,FE.nxU);
BC.USides        =  zeros(2,FE.nzU);
BC.WSides        =  zeros(2,FE.nzU);

BC.UTopBot(1,:)  =  U(1  ,:);
BC.UTopBot(2,:)  =  U(end,:);
BC.WTopBot(1,:)  =  W(1  ,:);
BC.WTopBot(2,:)  =  W(end,:);

BC.USides(1,:)   =  U(:,1  );
BC.USides(2,:)   =  U(:,end);
BC.WSides(1,:)   =  W(:,1  );
BC.WSides(2,:)   =  W(:,end);

if strcmp(CTX.MODEL.Type,'Erebus')
    x    =  (FE.CoordU(FE.MapUn(1,:),1)' - FE.W/2) ./ (CTX.MODEL.IOFlowWidth);
    ind  =  x >= -0.5 & x <= 0.5;
    if CTX.MODEL.IOFlowSym
        inflow  =  (cos(x(ind).*4*pi)+cos(x(ind).*2*pi))./2 .* -CTX.MODEL.IOFlowRate;
    else
        inflow  =  sin(x(ind).*2*pi) .* CTX.MODEL.IOFlowRate;
    end
    inflow             =  inflow - mean(inflow);
    inflow             =  inflow .* (1+cos(CTX.TIME.total.*2*pi/CTX.MODEL.IOFlowPeriod))./2;
    BC.WTopBot(2,ind)  =  BC.WTopBot(2,ind) + inflow;
    apply_botvz  =  1;
    apply_botvx  =  1;
end

%*****  optimized extracted boundary conditions  **************************

bc_ind   = [];
bc_val   = [];

if apply_topvx    ==  1
    
    ind      =  mapv(1,2:end-1,1)';
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nxU-2,1).*BC.UTopBot(1,2:end-1).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
    
end

if apply_topvz    ==  1
    
    ind      =  mapv(1,:,2)';
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nxU,1).*BC.WTopBot(1,:).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];

    
end

if apply_botvx    ==  1
    
    ind      =  mapv(end,2:end-1,1)';
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nxU-2,1).*BC.UTopBot(2,2:end-1).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
    
end

if apply_botvz    ==  1
    
    ind      =  mapv(end,:,2)';
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nxU,1).*BC.WTopBot(2,:).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
    
end

if apply_leftvx   ==  1
    
    ind      =  mapv(:,1,1);
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nzU,1).*BC.USides(1,:).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];

    
end

if apply_leftvz   ==  1
    
    ind      =  mapv(2:end-1,1,2);
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nzU-2,1).*BC.WSides(1,2:end-1).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
    
end

if apply_rightvx  ==  1
    
    ind      =  mapv(:,end,1);
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nzU,1).*BC.USides(2,:).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
    
end

if apply_rightvz  ==  1
    
    ind      =  mapv(2:end-1,end,2);
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nzU-2,1).*BC.WSides(2,2:end-1).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
    
end

% don't compute velocity solution if viscosity on upper cutoff
if strcmp(CTX.MODEL.Type,'Erebus')
    Mat = CTX.MP.Mat;
    Eta = CTX.MP.Eta;
    Max = 0.1.*CTX.RHEO.MaxEta;
    if strcmp(FE.ElType,'Q2Q1')
        indu     =  PQ1Q2(Eta,FE)>=Max & PQ1Q2(Mat,FE)<=1.05;
        indw     =  PQ1Q2(Eta,FE)>=Max & PQ1Q2(Mat,FE)<=1.05;
        indp     =        Eta    >=Max &       Mat    <=1.05;
    elseif strcmp(FE.ElType,'Q1Q1')
        indu     =        Eta    >=Max &       Mat    <=1.05;
        indw     =        Eta    >=Max &       Mat    <=1.05;
        indp     =        Eta    >=Max &       Mat    <=1.05;
    elseif strcmp(FE.ElType,'Q1P0')
        indu     =        Eta    >=Max &       Mat    <=1.05;
        indw     =        Eta    >=Max &       Mat    <=1.05;
        indp     =  PQ1El(Eta,FE)>=Max & PQ1El(Mat,FE)<=1.05;
    end
    bc_ind   =  [bc_ind;FE.DOFU(indu);FE.DOFW(indw);FE.DOFP(indp)];
    UU = U(:);  WW = W(:);
    tmpu     =  UU(indu)./diag(X(FE.DOFU(indu),FE.DOFU(indu)));
    tmpw     =  WW(indw)./diag(X(FE.DOFW(indw),FE.DOFW(indw)));
    tmpp     =  zeros(size(indp));
    bc_val   =  [bc_val;tmpu;tmpw;tmpp];
end

% fix P magnitude if no free surface
if apply_Pfix == 1
    P0       =  0;%(mean(CTX.MP.Rho(FE.MapP(:,round(FE.nxP/2)))) - CTX.SL.RhoRef).*CTX.PHYS.grav(2).*FE.D;
    ind      =  mapP(round(FE.nzP/2),round(FE.nxP/2));
    bc_ind   =  [bc_ind;ind(:)];
    tmp      =  P0.*ones(size(ind(:)))./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

[BC.ind,ind]     =  sort(bc_ind);
BC.val           =  bc_val(ind);

BC.free          =  1:length(L(1,:));
BC.free(BC.ind)  =  [];
TMP              =  L(:,BC.ind);
RHS              =  RHS - TMP*BC.val;
L                =  L(BC.free,BC.free);


end