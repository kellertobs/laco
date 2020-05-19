

function   Q2field  =  PQ1Q2(Q1field,FE)

n        =  length(Q1field(1,:));
Q2field  =  zeros(FE.NQ2,n);
[N,~]    =  ShapeFuncts([-1,0,1,-1,0,1,-1,0,1;-1,-1,-1,0,0,0,1,1,1],'Q1');


for k=1:n
    
    TMP                     =  Q1field(:,k);
    TMP                     =  TMP(FE.El2Q1)*N;
    Q2field(FE.El2Q2(:),k)  =  TMP(:);
    
end