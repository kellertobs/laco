

function   Q1field  =  PElQ1(Elfield,FE)

n        =  length(Elfield(1,:));
ElX      =  FE.CoordEl(:,1);
ElZ      =  FE.CoordEl(:,2);
Q1X      =  FE.CoordQ1(:,1);
Q1Z      =  FE.CoordQ1(:,2);
Q1field  =  zeros(FE.NQ1,n);
SZ       =  size(Q1field(:,1));

i   =  FE.El2Q1;
DX  =  abs(repmat(ElX,1,4) - Q1X(i));
DZ  =  abs(repmat(ElZ,1,4) - Q1Z(i));
DW  =  (1-DX./FE.hxQ1) .* (1-DZ./FE.hzQ1);

for k = 1:n
    
    El      =  repmat(Elfield(:,k),1,4);
    
    SUM     =  accumarray(i(:),El(:).*DW(:),SZ);
    WEIGHT  =  accumarray(i(:),       DW(:),SZ);
    
    Q1field(:,k)    =  SUM./WEIGHT;
    
end

end

