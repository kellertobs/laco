

function   Q1field  =  PQ2Q1(Q2field,FE)

n        =  length(Q2field(1,:));
Q1field  =  zeros(FE.NQ1,n);
[N,~]    =  ShapeFuncts([-1,1,-1,1;-1,-1,1,1],'Q2');


for k=1:n
    
    TMP                     =  Q2field(:,k);
    TMP                     =  TMP(FE.El2Q2)*N;
    Q1field(FE.El2Q1(:),k)  =  TMP(:);
    
end