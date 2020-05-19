

function   Elfield  =  PQ1El(Q1field,FE)

n        =  length(Q1field(1,:));
Elfield  =  zeros(FE.NEl,n);
[N,~]    =  ShapeFuncts([0;0],'Q1');


for k=1:n
    
    TMP           =  Q1field(:,k);
    TMP           =  TMP(FE.El2Q1)*N;
    Elfield(:,k)  =  TMP(:);
    
end