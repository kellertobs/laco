

function   IPfield  =  PQ2IP(Q2field,FE)

n        =  length(Q2field(1,:));
IPfield  =  zeros(FE.NIP,n);
N        =  FE.NiQ2;

for k=1:n
    
    TMP                     =  Q2field(:,k);
    TMP                     =  TMP(FE.El2Q2)*N;
    IPfield(FE.El2IP(:),k)  =  TMP(:);
    
end