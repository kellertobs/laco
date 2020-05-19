

function   Q2field  =  PElQ2(Elfield,FE)

n        =  length(Elfield(1,:));
ElX      =  FE.CoordEl(:,1);
ElZ      =  FE.CoordEl(:,2);
Q2X      =  FE.CoordQ2(:,1);
Q2Z      =  FE.CoordQ2(:,2);
Q2field  =  zeros(FE.NQ2,n);
SZ       =  size(Q2field(:,1));

i   =  FE.El2Q2;
DX  =  (repmat(ElX,1,9) - Q2X(i));
DZ  =  (repmat(ElZ,1,9) - Q2Z(i));
DW  =  max(1e-64,(1-DX./FE.hxQ2) .* (1-DZ./FE.hzQ2));

for k = 1:n
    
    El      =  repmat(Elfield(:,k),1,9);
    
    SUM     =  accumarray(i(:),El(:).*DW(:),SZ);
    WEIGHT  =  accumarray(i(:),       DW(:),SZ);
    
    Q2field(:,k)    =  SUM./WEIGHT;
    
end

end

