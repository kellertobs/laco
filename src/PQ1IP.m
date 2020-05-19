

function   IPfield  =  PQ1IP(Q1field,FE)

n        =  length(Q1field(1,:));
IPfield  =  zeros(FE.NIP,n);
N        =  FE.NiQ1;

for k=1:n
    
    TMP                     =  Q1field(:,k);
    TMP                     =  TMP(FE.El2Q1)*N;
    IPfield(FE.El2IP(:),k)  =  TMP(:);
    
end