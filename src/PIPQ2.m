
function  Q2field  =  PIPQ2(IPfield,FE)

n        =  length(IPfield(1,:));
Q2field  =  zeros(FE.NQ2,n);
N        =  FE.NiQ2;

for k=1:n
    
    TMP           =  IPfield(:,k);
    TMP           =  (TMP(FE.El2IP)*N') ./ (ones(size(FE.El2IP))*N');
    Q2field(:,k)  =  accumarray(FE.El2Q2(:),TMP(:)) ./ accumarray(FE.El2Q2(:),ones(size(TMP(:))));
    
end

end