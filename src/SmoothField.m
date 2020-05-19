% regularize field to avoid instability

function [field] = SmoothField(field,kappa,nsteps,FE,type)

ncomp = size(field,2);

switch type
    
    case 'El'
        
        map  =  FE.MapEl;
        nx   =  FE.nxEl;
        nz   =  FE.nzEl;
        
    case 'Q1'
        
        map  =  FE.MapQ1;
        nx   =  FE.nxQ1;
        nz   =  FE.nzQ1;
        
    case 'Q2'
        
        map  =  FE.MapQ2;
        nx   =  FE.nxQ2;
        nz   =  FE.nzQ2;
        
    case 'IP'
        
        map  =  FE.MapIP;
        nx   =  FE.nxIP;
        nz   =  FE.nzIP;
        
end

for n = 1:1:nsteps
    
    for k = 1:ncomp
        fk                   =  field(:,k);
        wk                   =  zeros(nz+2,nx+2);
        wk(2:end-1,2:end-1)  =  fk(map);
        wk(1  ,2:end-1)      =  fk(map(1  ,:));
        wk(end,2:end-1)      =  fk(map(end,:));
        wk(2:end-1,  1)      =  fk(map(:,  1));
        wk(2:end-1,end)      =  fk(map(:,end));
        
        diff        =  (wk(1:end-2,2:end-1) + wk(3:end,2:end-1) + wk(2:end-1,1:end-2) + wk(2:end-1,3:end)) - 4.*wk(2:end-1,2:end-1);
        fk(map)     =  fk(map) + kappa/8 .* diff;
        field(:,k)  =  fk;
    end
    
end

end





