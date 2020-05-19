% Function returns coordinates and weights for numerical integration
% points in natural coordinate system (local element)
%
% [FE]  =  IntCoords(FE)
% Function returns local coordinates of integration points PAR.Si (2,NIP)
% and integration weights PAR.Wi (1,NIP);
%
% sorting of integration points is
%
%                   4
% NIP  =  9:  1 o---o---o 7
%               |   |5  |
%             2 o---o---o 8
%               |   |   |
%             3 o---o---o 9
%                   6   
%
% created  20140729 Tobias Keller
% modified 20190418 Tobias Keller


function  [Si,Wi]  =  IntCoords(FE)

%***  prepare arrays

Si      =  zeros(2,FE.IPpEl);
Wi      =  zeros(1,FE.IPpEl);


%***  local coordinates of integration points for NIP = 9

c        =  sqrt(3/5);

Si(:,1)  =  [-c; -c];
Si(:,2)  =  [ 0; -c];
Si(:,3)  =  [ c; -c];

Si(:,4)  =  [-c;  0];
Si(:,5)  =  [ 0;  0];
Si(:,6)  =  [ c;  0];

Si(:,7)  =  [-c;  c];
Si(:,8)  =  [ 0;  c];
Si(:,9)  =  [ c;  c];


%***  integration weights of integration points for NIP = 9

Wi(:,1)  =  5/9*5/9;
Wi(:,2)  =  5/9*8/9;
Wi(:,3)  =  5/9*5/9;

Wi(:,4)  =  8/9*5/9;
Wi(:,5)  =  8/9*8/9;
Wi(:,6)  =  8/9*5/9;

Wi(:,7)  =  5/9*5/9;
Wi(:,8)  =  5/9*8/9;
Wi(:,9)  =  5/9*5/9;


end


