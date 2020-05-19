% PrintDiagnostics    EDIFICE: Report model diagnostics
%
% []  =  PrintDiagnostics(CTX)
%
%   Function reports model diagnostics by printing to standard output
%
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function  []  =  PrintDiagnostics(CTX)

fprintf('    max eta      = ')
fprintf(1,'%e',max(CTX.MP.EtaVEP));
fprintf('    min eta      = ')
fprintf(1,'%e\n',min(CTX.MP.EtaVEP));

fprintf('    max stress   = ')
fprintf(1,'%e',max(CTX.MP.TII(:,1)));
fprintf('    min stress   = ')
fprintf(1,'%e\n',min(CTX.MP.TII(:,1)));

fprintf('    max strainr  = ')
fprintf(1,'%e',max(CTX.MP.EII(:,1)));
fprintf('    min strainr  = ')
fprintf(1,'%e\n',min(CTX.MP.EII(:,1)));

fprintf('    max U        = ')
fprintf(1,'%e',max(abs(CTX.SL.U))*CTX.TIME.spyr);
fprintf('    max W        = ')
fprintf(1,'%e\n',max(abs(CTX.SL.W))*CTX.TIME.spyr);

fprintf('    max deform   = ')
fprintf(1,'%f',CTX.FE.MaxDef );
fprintf('    min deform   = ')
fprintf(1,'%f\n',CTX.FE.MinDef );
fprintf('    max volume   = ')
fprintf(1,'%f',max(CTX.FE.ElVol)/mean(CTX.FE.ElVol) );
fprintf('    min volume   = ')
fprintf(1,'%f\n\n',min(CTX.FE.ElVol)/mean(CTX.FE.ElVol) );

end