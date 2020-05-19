

function   Elfield  =  PQ2El(Q2field,FE)

n        =  length(Q2field(1,:));
Elfield  =  zeros(FE.NEl,n);
[N,~]    =  ShapeFuncts([0;0],'Q2');


for k=1:n
    
    TMP           =  Q2field(:,k);
    TMP           =  TMP(FE.El2Q2)*N;
    Elfield(:,k)  =  TMP(:);
    
end