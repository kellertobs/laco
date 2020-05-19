% AssembleOperator    EDIFICE: Construct coefficient matrix
%
% [L,RHS]  =  AssembleOperator(MP,SL,CTX)
%
%   Function constructs coefficient matrix L and right-hand side vector R
%   for linear Stokes problem according to problem dimensions and shape 
%   functions in FE, and material point properties in MP.
%
%   created 20161115 Tobias Keller
%   modified  20170508  Tobias Keller
%   modified  20200227   Tobias Keller


function [L,RHS]  =  AssembleOperator(MP,SL,CTX)

FE  =  CTX.FE;


%*****  prepare some parameters  ******************************************

g       =  CTX.PHYS.grav;
EtaVP   =  MP.EtaVP  (FE.El2IP);
EtaVEP  =  MP.EtaVEP (FE.El2IP);
Rho     =  MP.Rho    (FE.El2IP);
Src     =  MP.SrcdVdt(FE.El2IP);

if strcmp(CTX.RHEO.Elasticity,'ON')
    Txxo     =  MP.Taur(:,1);
    Tzzo     =  MP.Taur(:,2);
    Txzo     =  MP.Taur(:,3);
    ChiTxxo  =  MP.Chi(FE.El2IP).*Txxo(FE.El2IP);
    ChiTzzo  =  MP.Chi(FE.El2IP).*Tzzo(FE.El2IP);
    ChiTxzo  =  MP.Chi(FE.El2IP).*Txzo(FE.El2IP);
end

NiU      =  FE.NiU;
NiP      =  FE.NiP;
dNdSiU   =  FE.dNdSiU;
dNdSiP   =  FE.dNdSiP;
Wi       =  FE.Wi;

el2u     =  FE.El2Un;
el2p     =  FE.El2Pn;
el2v     =  FE.El2Vd;

NEl      =  FE.NEl;
nu       =  FE.UpEl;
np       =  FE.PpEl;
nip      =  FE.IPpEl;
nv       =  nu*2;

nelblo   =  FE.nblock;
nelblo   =  min(NEl, nelblo);
nblo     =  ceil(NEl/nelblo);


%*****  allocate arrays  **************************************************

VV_all    =  zeros(NEl,nv*nv);
VP_all    =  zeros(NEl,nv*np);
PP_all    =  zeros(NEl,np*np);
RhsV_all  =  zeros(NEl,nv);
RhsP_all  =  zeros(NEl,np);


%*****  block assembly loop  **********************************************

il        =  1;
iu        =  nelblo;

for ib = 1:nblo
    
    ind     =  il:iu;
    
    COORDX  =  reshape(FE.CoordU(el2u(ind,:),1),nelblo,nu);
    COORDZ  =  reshape(FE.CoordU(el2u(ind,:),2),nelblo,nu);
    
    
    %*****  allocate arrays  **************************************
    
    VV_block    =  zeros(nelblo,nv*nv);
    VP_block    =  zeros(nelblo,nv*np);
    PP_block    =  zeros(nelblo,np*np);
    RhsV_block  =  zeros(nelblo,nv);
    RhsP_block  =  zeros(nelblo,np);
    invJx       =  zeros(nelblo,2);
    invJz       =  zeros(nelblo,2);

    
    %*****  integration loop  *************************************
    
    for ip=1:nip
        
        NU            =  NiU(:,ip);
        NP            =  NiP(:,ip);
        dNdSU         =  dNdSiU(:,:,ip);
        dNdSP         =  dNdSiP(:,:,ip);
        
        Jx            =  COORDX*dNdSU';
        Jz            =  COORDZ*dNdSU';
        detJ          =  Jx(:,1).*Jz(:,2) - Jx(:,2).*Jz(:,1);
        
        invdetJ       =  1.0./detJ;
        invJx(:,1)    =  +Jz(:,2).*invdetJ;
        invJx(:,2)    =  -Jz(:,1).*invdetJ;
        invJz(:,1)    =  -Jx(:,2).*invdetJ;
        invJz(:,2)    =  +Jx(:,1).*invdetJ;
        
        dNdxU         =  invJx*dNdSU;
        dNdzU         =  invJz*dNdSU;
        
        dNdxP         =  invJx*dNdSP;
        dNdzP         =  invJz*dNdSP;
        
        NNU           =  ones(nelblo,1)*NU.';
        NNP           =  ones(nelblo,1)*NP.';
        W             =  Wi(ip).*detJ;
        WU            =  W(:,ones(1,nu));
        WP            =  W(:,ones(1,np));
        
        
        %*****  VV matrix  ****************************************
        
        a0    =  2/3.*EtaVEP(ind,ip);
        a1    = -1/3.*EtaVEP(ind,ip);
        a2    =  1/2.*EtaVEP(ind,ip);
        
        indx  = 1;
        for i = 1:nu
            % x-velocity equation
            for j = 1:nu
                VV_block(:,indx) = VV_block(:,indx) - (a0.*dNdxU(:,i).*dNdxU(:,j) ...
                                                    +  a2.*dNdzU(:,i).*dNdzU(:,j)).*W;
                indx = indx+1;
                VV_block(:,indx) = VV_block(:,indx) - (a1.*dNdxU(:,i).*dNdzU(:,j) ...
                                                    +  a2.*dNdzU(:,i).*dNdxU(:,j)).*W;
                indx = indx+1;
            end
            % y-velocity equation
            for j = 1:nu
                VV_block(:,indx) = VV_block(:,indx) - (a1.*dNdzU(:,i).*dNdxU(:,j) ...
                                                    +  a2.*dNdxU(:,i).*dNdzU(:,j)).*W;
                indx = indx+1;
                VV_block(:,indx) = VV_block(:,indx) - (a0.*dNdzU(:,i).*dNdzU(:,j) ...
                                                    +  a2.*dNdxU(:,i).*dNdxU(:,j)).*W;
                indx = indx+1;
            end
        end
        
        
        %***** SL matrix  *****************************************
        
        for i=1:np
            VP_block(:,(i-1)*nv+(1:2:nv)) =  VP_block(:,(i-1)*nv+(1:2:nv)) + dNdxU.*NNP(:,i).*WU;
            VP_block(:,(i-1)*nv+(2:2:nv)) =  VP_block(:,(i-1)*nv+(2:2:nv)) + dNdzU.*NNP(:,i).*WU;
        end
        
        
        %*****  RhsV vector  **************************************
        
        RhsV_block(:,1:2:nv) = RhsV_block(:,1:2:nv) - (0.*(Rho(ind,ip)-SL.RhoRef)).*NNU.*WU;
        RhsV_block(:,2:2:nv) = RhsV_block(:,2:2:nv) - (g.*(Rho(ind,ip)-SL.RhoRef)).*NNU.*WU;
        
        if strcmp(CTX.RHEO.Elasticity,'ON')
            RhsV_block(:,1:2:nv) =  RhsV_block(:,1:2:nv) + (dNdxU.*ChiTxxo(ind,ip) + dNdzU.*ChiTxzo(ind,ip)).*WU;
            RhsV_block(:,2:2:nv) =  RhsV_block(:,2:2:nv) + (dNdxU.*ChiTxzo(ind,ip) + dNdzU.*ChiTzzo(ind,ip)).*WU;
        end
        
        
        %*****  PP matrix / RhsP vector ************************
        
%         stab  =  CTX.SL.StabFact.*FE.ElVol(ind)./EtaVEP(ind,ip);
        stab  =  CTX.SL.StabFact;

        indx   =  1;
        for i = 1:np
            for j = 1:np
                PP_block(:,indx)  =  PP_block(:,indx) + (dNdxP(:,i).*dNdxP(:,j)  ...
                                                      +  dNdzP(:,i).*dNdzP(:,j)) ... 
                                                      .* stab .* W;
%                 PP_block(:,indx)  =  PP_block(:,indx) + stab .* NNP .* WP;
                indx = indx+1;
            end
        end
        
        RhsP_block  =  RhsP_block + Src(ind,ip) .* NNP .* WP;
        
    end
    
    
    %*****  collect element matrices  *****************************
    
    VV_all(ind,:)    =  VV_block;
    VP_all(ind,:)    =  VP_block;
    PP_all(ind,:)    =  PP_block;
    RhsV_all(ind,:)  =  RhsV_block;
    RhsP_all(ind,:)  =  RhsP_block;

    
    %*****  adjust blocksize and arrays  **********************************
    
    il  =  il+nelblo;
    if(ib==nblo-1)
        nelblo = NEl-iu;
    end
    
    iu  =  iu+nelblo;
    
end


%*****  create matrix indices  ********************************************

%***  VV matrix
indx_j  =  repmat(1:nv,nv,1); indx_i = indx_j';
VV_i    =  el2v(:,indx_i(:));
VV_j    =  el2v(:,indx_j(:));

%***  PP matrix
indx_j  =  repmat(1:np,np,1); indx_i = indx_j';
PP_i    =  el2p(:,indx_i(:));
PP_j    =  el2p(:,indx_j(:));

%***  VP matrix
indx_i  =  repmat(1:nv,np,1)';
indx_j  =  repmat(1:np,nv,1);
VP_i    =  el2v(:,indx_i(:));
VP_j    =  el2p(:,indx_j(:));

%*****  assemble global block matrices  ***********************************

VV     =  sparse(VV_i(:),  VV_j(:),  VV_all(:));
PP     =  sparse(PP_i(:),  PP_j(:),  PP_all(:));
VP     =  sparse(VP_i(:),  VP_j(:),  VP_all(:));
PV     =  VP.';

%*****  assemble RHS block vectors  ***************************************

RhsV  =  accumarray(el2v(:), RhsV_all(:));
RhsP  =  accumarray(el2p(:), RhsP_all(:));


%***  assemble global operator  *******************************************

L     =  [VV   VP ; ...
          PV   PP ];

RHS   =  [RhsV; RhsP];


end


