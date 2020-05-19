% NonlinearSolver    EDIFICE: Calculates solution to non-linear Stokes equations
%
% [CTX] = NonlinearSolver(MP,CTX)
%
%   Function finds solution to non-linear governing equations. It operates
%   according to options set in CTX.SL; it iteratively calls LinearSolver()
%   to update the solution given a set of non-linear coefficients, and then
%   updates these to the new solution guess in UpdateMaterialPoints().
%
%   created 20161115 Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20170508  Tobias Keller
%   modified  20200227   Tobias Keller


function  [CTX] = NonlinearSolver(CTX)
    
CTX.SLo = CTX.SL;
CTX.MPo = CTX.MP;
CTX.FEo = CTX.FE;

CTX.SL.it     =  0;
CTX.SL.fnorm  =  1;
fnorm0        =  1;


while (  CTX.SL.fnorm        > CTX.SL.atol   ... 
      && CTX.SL.fnorm/fnorm0 > CTX.SL.rtol   ...
      && CTX.SL.it           < CTX.SL.maxits )
    
    CTX.SL.it = CTX.SL.it+1;
    
    if CTX.SL.it < 10
        fprintf(1,'    %i',CTX.SL.it)
    else
        fprintf(1,'   %i',CTX.SL.it)
    end
    
    if strcmp(CTX.FE.LagrMesh,'ON')
        CTX  =  AdvectFE(CTX);
    end
    
    CTX  =  LinearSolver(CTX);
    CTX  =  UpdateMaterialPoints(CTX);
    
    if CTX.SL.it == 1 || CTX.SL.fnorm > fnorm0
        fnorm0 = CTX.SL.fnorm;
    end
    
end

