% LivePlotting    EDIFICE: Save model results to file
%
% []  =  SaveToFile(CTX)
%
%   Function saves model solution and material point properties to file.
%
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function  []  =  SaveToFile(CTX)

%*****  save CTX to file  *****************************************
fprintf(1,'******************************************************\n')
fprintf(1,'***************     save frame # %d     **************\n',CTX.IO.frame)
fprintf(1,'******************************************************\n\n')

CTX = rmfield(CTX,'SLo');
CTX = rmfield(CTX,'MPo');
CTX = rmfield(CTX,'FEo');
CTX.SL = rmfield(CTX.SL,'S' );

ContName  =  [CTX.IO.DataDir '/' CTX.IO.RunID '/' CTX.IO.RunID '_cont.mat'];
SaveName  =  [CTX.IO.DataDir '/' CTX.IO.RunID '/' CTX.IO.RunID '_' num2str(CTX.IO.frame) '.mat'];

save(ContName,'CTX')
save(SaveName,'CTX')

end