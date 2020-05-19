
function  Q1field  =  PIPQ1(IPfield,FE)

n        =  length(IPfield(1,:));
Q1field  =  zeros(FE.NQ1,n);
N        =  FE.NiQ1;

for k=1:n
    
    TMP           =  IPfield(:,k);
    TMP           =  (TMP(FE.El2IP)*N') ./ (ones(size(FE.El2IP))*N');
    Q1field(:,k)  =  accumarray(FE.El2Q1(:),TMP(:)) ./ accumarray(FE.El2Q1(:),ones(size(TMP(:))));
    
end

end