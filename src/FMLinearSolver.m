
function  [CTX] = LinearSolver(CTX)

SL  =  CTX.SL;
MP  =  CTX.MP;
BC  =  CTX.BC;
FE  =  CTX.FE;
S   =  SL.SFM;
X   =  SL.XFM;


%*****  assemble operator and right hand side vector  *********************

[L,RHS]     =  FMAssembleOperator(MP,SL,CTX);

if CTX.SL.outer_it<=1
    X       =  max(1e-16,min(1e16,sqrt(abs(diag(L)))));
    X       =  diag(sparse(1./X));
    
    L       =  X*L*X;
    RHS     =  X*RHS;
else
    L    =  X*L*X;
    RHS  =  X*RHS;
end

[L,RHS,BC]  =  FMApplyBC(L,RHS,X,BC,FE,CTX);


%*****  get non-linear residual  ******************************************

F           =  zeros(size(S));
F(BC.free)  =  L*(1./diag(X(BC.free,BC.free)).*S(BC.free)) - RHS(BC.free);
F(BC.ind)   =  0.;

SL.F.U      =  F(FE.DOFU);
SL.F.W      =  F(FE.DOFW);
SL.F.P      =  F(FE.DOFP);

SL.FMfnorm  =  norm(F,2) / max(norm(RHS,2),1.e-32);

if SL.outer_it == 1 || SL.FMfnorm > SL.FMfnorm0
    SL.FMfnorm0  =  SL.FMfnorm;
end

fprintf(1,'   FM abs =')
fprintf(1,' %4.4e',SL.FMfnorm)
fprintf(1,'  rel =')
fprintf(1,' %4.4e\n',SL.FMfnorm./SL.FMfnorm0)


if SL.FMfnorm < SL.outer_atol && SL.outer_it > 1
    CTX.SL  =  SL;
    CTX.BC  =  BC;
    return
end

if SL.FMfnorm > 1.e+20 && SL.outer_it > 1
    SaveToFile(CTX);
    error('!!!  Error: Diverging solution, stop and try again  !!!')
end


%*****  solve linear system of equations  *********************************

S           =  zeros(size(S));
S(BC.free)  =  L \ RHS(BC.free);
S(BC.ind)   =  BC.val;

SL.SFM      =  X*S;

SL.U  =  SL.SFM(FE.DOFU);
SL.W  =  SL.SFM(FE.DOFW);
SL.P  =  SL.SFM(FE.DOFP);

SL.Plith  =  CTX.SL.RhoRef.*CTX.PHYS.grav(:,2).*(FE.CoordP(:,2)-min(FE.CoordP(CTX.FE.MapPn(1,:),2))) + CTX.BC.SurfPres;
SL.Pt     =  SL.P  + SL.Plith;

CTX.SL  =  SL;
CTX.BC  =  BC;

end

