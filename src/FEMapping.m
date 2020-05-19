% FeMapping    Create mapping arrays for use in EREBUS
%
% [FE] = FeMapping(FE)
%
%   creates mapping arrays mapping between elements and degrees of freedom 
%   on the FE mesh
%
%   created   20140729  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20190418  Tobias Keller


function [FE] = FEMapping(FE)


%*****  create 2D mapping arrays  *****************************************

mapEl  =  reshape(1:FE.NEl,FE.nzEl,FE.nxEl);
mapQ1  =  reshape(1:FE.NQ1,FE.nzQ1,FE.nxQ1);
mapQ2  =  reshape(1:FE.NQ2,FE.nzQ2,FE.nxQ2);
mapUn  =  reshape(1:FE.NU ,FE.nzU ,FE.nxU );
mapPn  =  reshape(1:FE.NP ,FE.nzP ,FE.nxP );
mapU   =  reshape(1:2:2*FE.NU,FE.nzU,FE.nxU);
mapW   =  reshape(2:2:2*FE.NU,FE.nzU,FE.nxU);
mapV   =  reshape([mapU,mapW],FE.nzU,FE.nxU,2);
mapP   =  reshape(1:FE.NP,FE.nzP,FE.nxP) + 2*FE.NU;
mapIP  =  reshape(1:FE.NIP,FE.nzIP,FE.nxIP);


%*****  compose El-to-DOF/IP mapping arrays  ******************************

El2P0     =  zeros(FE.NEl,  1       );
El2Q1     =  zeros(FE.NEl,  FE.Q1pEl);
El2Q2     =  zeros(FE.NEl,  FE.Q2pEl);
if FE.NU == FE.NQ1
    El2V      =  zeros(FE.NEl,2*FE.Q1pEl);
else
    El2V      =  zeros(FE.NEl,2*FE.Q2pEl);
end
El2IP     =  zeros(FE.NEl,  FE.IPpEl);

iel = 1;
for i = 1:1:FE.nxEl
    for j = 1:1:FE.nzEl

        El2P0(iel,:)   = mapEl(iel);

        iq1  =  i;   jq1  =  j;           
        El2Q1(iel,:)   = [mapQ1(jq1,iq1  ),mapQ1(jq1+1,iq1  ), ...
                          mapQ1(jq1,iq1+1),mapQ1(jq1+1,iq1+1)  ];
                      
        iq2  =  2*i-1;   jq2  =  2*j-1;           
        El2Q2(iel,:)   = [mapQ2(jq2,iq2  ),mapQ2(jq2+1,iq2  ), mapQ2(jq2+2,iq2  ), ...
                          mapQ2(jq2,iq2+1),mapQ2(jq2+1,iq2+1), mapQ2(jq2+2,iq2+1), ...                  
                          mapQ2(jq2,iq2+2),mapQ2(jq2+1,iq2+2), mapQ2(jq2+2,iq2+2)  ];
        
        if FE.NU == FE.NQ1
            iv = iq1; jv = jq1;
            El2V(iel,:)    = [mapV(jv,iv  ,1),mapV(jv,iv,  2),mapV(jv+1,iv  ,1),mapV(jv+1,iv  ,2), ...
                              mapV(jv,iv+1,1),mapV(jv,iv+1,2),mapV(jv+1,iv+1,1),mapV(jv+1,iv+1,2)  ];
        else
            iv = iq2; jv = jq2;
            El2V(iel,:)    = [mapV(jv,iv  ,1),mapV(jv,iv,  2),mapV(jv+1,iv  ,1),mapV(jv+1,iv  ,2),mapV(jv+2,iv  ,1),mapV(jv+2,iv  ,2), ...
                              mapV(jv,iv+1,1),mapV(jv,iv+1,2),mapV(jv+1,iv+1,1),mapV(jv+1,iv+1,2),mapV(jv+2,iv+1,1),mapV(jv+2,iv+1,2), ...
                              mapV(jv,iv+2,1),mapV(jv,iv+2,2),mapV(jv+1,iv+2,1),mapV(jv+1,iv+2,2),mapV(jv+2,iv+2,1),mapV(jv+2,iv+2,2)  ];
        end

        iip =  3*i-2;   jip =  3*j-2;
        El2IP(iel,:)   = [mapIP(jip,iip  ),mapIP(jip+1,iip  ), mapIP(jip+2,iip  ), ...
                          mapIP(jip,iip+1),mapIP(jip+1,iip+1), mapIP(jip+2,iip+1), ...                  
                          mapIP(jip,iip+2),mapIP(jip+1,iip+2), mapIP(jip+2,iip+2)  ];

        iel = iel + 1;
    end
end
    

%***  return structures
FE.MapEl  =  mapEl;
FE.MapQ1  =  mapQ1;
FE.MapQ2  =  mapQ2;
FE.MapIP  =  mapIP;
FE.MapUn  =  mapUn;
FE.MapPn  =  mapPn;
FE.MapU   =  mapU;
FE.MapW   =  mapW;
FE.MapV   =  mapV;
FE.MapP   =  mapP;
FE.MapT   =  mapQ1;
FE.El2Q1  =  El2Q1;
FE.El2Q2  =  El2Q2;
FE.El2IP  =  El2IP;
FE.El2Vd  =  El2V;
FE.El2Ud  =  El2V(:,1:2:end);
FE.El2Wd  =  El2V(:,2:2:end);
switch FE.ElType
    case 'Q2Q1'
        FE.El2Un  =  El2Q2;
        FE.El2Pn  =  El2Q1;
        FE.El2Pd  =  El2Q1 + 2*FE.NQ2;
    case 'Q1Q1'
        FE.El2Un  =  El2Q1;
        FE.El2Pn  =  El2Q1;
        FE.El2Pd  =  El2Q1 + 2*FE.NQ1;
    case 'Q1P0'
        FE.El2Un  =  El2Q1;
        FE.El2Pn  =  El2P0;
        FE.El2Pd  =  El2P0 + 2*FE.NQ1;
end
FE.El2T   =  El2Q1;


