% regularize field to avoid instability

function [field] = ExplicitAdvDiffSolver(field,oldfield,vel,oldvel,dv,olddv,src,srco,kappa,dt,type,CTX)

FE     =  CTX.FE;
ncomp  =  size(field,2);

switch type
    
    case 'El'
        
        h       =  FE.hzQ1;
        map     =  FE.MapEl;
        nx      =  FE.nxEl;
        nz      =  FE.nzEl;
        if FE.NU == FE.NQ2
            vel     =  PQ2El(vel,FE);
            oldvel  =  PQ2El(oldvel,FE);
        else
            vel     =  PQ1El(vel,FE);
            oldvel  =  PQ1El(oldvel,FE);
        end

    case 'Q1'
        
        h    =  FE.hzQ1;
        map  =  FE.MapQ1;
        nx   =  FE.nxQ1;
        nz   =  FE.nzQ1;
        if FE.NU == FE.NQ2
            vel     =  PQ2Q1(vel,FE);
            oldvel  =  PQ2Q1(oldvel,FE);
        end
        
    case 'Q2'
        
        h    =  FE.hzQ2;
        map  =  FE.MapQ2;
        nx   =  FE.nxQ2;
        nz   =  FE.nzQ2;
        if FE.NU == FE.NQ1
            vel     =  PQ1Q2(vel,FE);
            oldvel  =  PQ1Q2(oldvel,FE);
        end
        
end

ic   =  3:nx+2;
im   =  2:nx+1;
ip   =  4:nx+3;
imm  =  1:nx+0;
ipp  =  5:nx+4;

jc   =  3:nz+2;
jm   =  2:nz+1;
jp   =  4:nz+3;
jmm  =  1:nz+0;
jpp  =  5:nz+4;

for k = 1:ncomp
    
    fieldk   =     field(:,k);
    fieldko  =  oldfield(:,k);
    
    a        =  ghost_field(fieldk ,nx,nz,map);
    ao       =  ghost_field(fieldko,nx,nz,map);

    lapl     =  -(4.*a (jc,ic) - a (jc,im) - a (jc,ip) - a (jm,ic) - a (jp,ic))./h^2;
    laplo    =  -(4.*ao(jc,ic) - ao(jc,im) - ao(jc,ip) - ao(jm,ic) - ao(jp,ic))./h^2;
    
    if strcmp(CTX.SL.Advection(1:3),'UPW')
        
        velx     =  vel(:,1);
        velz     =  vel(:,2);
        dvx      =  dv(:,1);
        dvz      =  dv(:,2);
        
        u        =  ghost_field(velx,nx,nz,map);
        w        =  ghost_field(velz,nx,nz,map);
        du       =  ghost_field(dvx ,nx,nz,map);
        dw       =  ghost_field(dvz ,nx,nz,map);
        
        u = u + du;  w = w + dw;
        
        um       =  min(u(jc,ic),0);
        up       =  max(u(jc,ic),0);
        wm       =  min(w(jc,ic),0);
        wp       =  max(w(jc,ic),0);
        
        velxo    =  oldvel(:,1);
        velzo    =  oldvel(:,2);
        dvxo     =  olddv(:,1);
        dvzo     =  olddv(:,2);
        
        uo       =  ghost_field(velxo,nx,nz,map);
        wo       =  ghost_field(velzo,nx,nz,map);
        duo      =  ghost_field(dvxo ,nx,nz,map);
        dwo      =  ghost_field(dvzo ,nx,nz,map);
        
        uo = uo + duo;  wo = wo + dwo;
        
        umo      =  min(uo(jc,ic),0);
        upo      =  max(uo(jc,ic),0);
        wmo      =  min(wo(jc,ic),0);
        wpo      =  max(wo(jc,ic),0);
        
        if strcmp(CTX.SL.Advection(4),'2')
            
            da_dzm = ( 3 * a(jc,ic) - 4 * a(jm,ic) + a(jmm,ic))./(2 * h);
            da_dzp = (-3 * a(jc,ic) + 4 * a(jp,ic) - a(jpp,ic))./(2 * h);
            da_dxm = ( 3 * a(jc,ic) - 4 * a(jc,im) + a(jc,imm))./(2 * h);
            da_dxp = (-3 * a(jc,ic) + 4 * a(jc,ip) - a(jc,ipp))./(2 * h);
            
            adv     =  um .* da_dxp + up .* da_dxm + wm .* da_dzp + wp .* da_dzm;
            
            da_dzm = ( 3 * ao(jc,ic) - 4 * ao(jm,ic) + ao(jmm,ic))./(2 * h);
            da_dzp = (-3 * ao(jc,ic) + 4 * ao(jp,ic) - ao(jpp,ic))./(2 * h);
            da_dxm = ( 3 * ao(jc,ic) - 4 * ao(jc,im) + ao(jc,imm))./(2 * h);
            da_dxp = (-3 * ao(jc,ic) + 4 * ao(jc,ip) - ao(jc,ipp))./(2 * h);
            
            advo    =  umo .* da_dxp + upo .* da_dxm + wmo .* da_dzp + wpo .* da_dzm;
            
        elseif strcmp(CTX.SL.Advection(4),'3')
            
            da_dxm  =  ( 2 * a(jc,im) + 3 * a(jc,ic) - 6 * a(jc,im) + a(jc,imm))./(6 * h);
            da_dxp  =  (-2 * a(jc,ip) - 3 * a(jc,ic) + 6 * a(jc,ip) - a(jc,ipp))./(6 * h);
            da_dzm  =  ( 2 * a(jp,ic) + 3 * a(jc,ic) - 6 * a(jm,ic) + a(jmm,ic))./(6 * h);
            da_dzp  =  (-2 * a(jm,ic) - 3 * a(jc,ic) + 6 * a(jp,ic) - a(jpp,ic))./(6 * h);
            
            adv     =  um .* da_dxp + up .* da_dxm + wm .* da_dzp + wp .* da_dzm;
            
            da_dxm  =  ( 2 * ao(jc,im) + 3 * ao(jc,ic) - 6 * ao(jc,im) + ao(jc,imm))./(6 * h);
            da_dxp  =  (-2 * ao(jc,ip) - 3 * ao(jc,ic) + 6 * ao(jc,ip) - ao(jc,ipp))./(6 * h);
            da_dzm  =  ( 2 * ao(jp,ic) + 3 * ao(jc,ic) - 6 * ao(jm,ic) + ao(jmm,ic))./(6 * h);
            da_dzp  =  (-2 * ao(jm,ic) - 3 * ao(jc,ic) + 6 * ao(jp,ic) - ao(jpp,ic))./(6 * h);
            
            advo    =  umo .* da_dxp + upo .* da_dxm + wmo .* da_dzp + wpo .* da_dzm;
            
        end
        
        divdv  = (dw(jp,ic)-dw(jm,ic))./h + (du(jc,ip)-du(jc,im))./h;
        adv    = adv + divdv .* a(jc,ic);
          
        divdvo = (dwo(jp,ic)-dwo(jm,ic))./h + (duo(jc,ip)-duo(jc,im))./h;
        advo   = advo + divdvo .* ao(jc,ic);
        
    elseif strcmp(CTX.SL.Advection(1:3),'FRM')
        
        velx     =  vel(:,1);
        velz     =  vel(:,2);
        dvx      =  dv(:,1);
        dvz      =  dv(:,2);
        
        u        =  ghost_field(velx,nx,nz,map);
        w        =  ghost_field(velz,nx,nz,map);
        du       =  ghost_field(dvx ,nx,nz,map);
        dw       =  ghost_field(dvz ,nx,nz,map);
        
        velxo    =  oldvel(:,1);
        velzo    =  oldvel(:,2);
        dvxo     =  olddv(:,1);
        dvzo     =  olddv(:,2);
        
        uo       =  ghost_field(velxo,nx,nz,map);
        wo       =  ghost_field(velzo,nx,nz,map);
        duo      =  ghost_field(dvxo ,nx,nz,map);
        dwo      =  ghost_field(dvzo ,nx,nz,map);
        
        wp = (w(jp,ic)+w(jc,ic))./2; 
        wm = (w(jm,ic)+w(jc,ic))./2;
        up = (u(jc,ip)+u(jc,ic))./2;  
        um = (u(jc,im)+u(jc,ic))./2;
        divv = (wp-wm)./h + (up-um)./h;
        
        wp = (w(jp,ic)+dw(jp,ic)+w(jc,ic)+dw(jc,ic))./2;
        wm = (w(jm,ic)+dw(jm,ic)+w(jc,ic)+dw(jc,ic))./2;
        up = (u(jc,ip)+du(jc,ip)+u(jc,ic)+du(jc,ic))./2;  
        um = (u(jc,im)+du(jc,im)+u(jc,ic)+du(jc,ic))./2;
        
        acc = a(jc,ic);
        ajp = a(jp,ic); ajpp = a(jpp,ic); 
        ajm = a(jm,ic); ajmm = a(jmm,ic);
        aip = a(jc,ip); aipp = a(jc,ipp); 
        aim = a(jc,im); aimm = a(jc,imm);
        
        adv  =  ((up .*(-aipp + 5.*(aip+acc)-aim )./8 - abs(up).*(-aipp + 3.*(aip-acc)+aim )./8) - ...
                 (um .*(-aip  + 5.*(acc+aim)-aimm)./8 - abs(um).*(-aip  + 3.*(acc-aim)+aimm)./8))./h ...
              + ((wp .*(-ajpp + 5.*(ajp+acc)-ajm )./8 - abs(wp).*(-ajpp + 3.*(ajp-acc)+ajm )./8) - ...
                 (wm .*(-ajp  + 5.*(acc+ajm)-ajmm)./8 - abs(wm).*(-ajp  + 3.*(acc-ajm)+ajmm)./8))./h;
        adv  =  adv - acc.*divv;
             
        wp = (wo(jp,ic)+wo(jc,ic))./2; 
        wm = (wo(jm,ic)+wo(jc,ic))./2;
        up = (uo(jc,ip)+uo(jc,ic))./2;  
        um = (uo(jc,im)+uo(jc,ic))./2;
        divv = (wp-wm)./h + (up-um)./h;
        
        wp = (wo(jp,ic)+dwo(jp,ic)+wo(jc,ic)+dwo(jc,ic))./2;
        wm = (wo(jm,ic)+dwo(jm,ic)+wo(jc,ic)+dwo(jc,ic))./2;
        up = (uo(jc,ip)+duo(jc,ip)+uo(jc,ic)+duo(jc,ic))./2;  
        um = (uo(jc,im)+duo(jc,im)+uo(jc,ic)+duo(jc,ic))./2;
        
        acc = ao(jc,ic);
        ajp = ao(jp,ic); ajpp = ao(jpp,ic); 
        ajm = ao(jm,ic); ajmm = ao(jmm,ic);
        aip = ao(jc,ip); aipp = ao(jc,ipp); 
        aim = ao(jc,im); aimm = ao(jc,imm);
        
        advo =  ((up .*(-aipp + 5.*(aip+acc)-aim )./8 - abs(up).*(-aipp + 3.*(aip-acc)+aim )./8) - ...
                 (um .*(-aip  + 5.*(acc+aim)-aimm)./8 - abs(um).*(-aip  + 3.*(acc-aim)+aimm)./8))./h ...
              + ((wp .*(-ajpp + 5.*(ajp+acc)-ajm )./8 - abs(wp).*(-ajpp + 3.*(ajp-acc)+ajm )./8) - ...
                 (wm .*(-ajp  + 5.*(acc+ajm)-ajmm)./8 - abs(wm).*(-ajp  + 3.*(acc-ajm)+ajmm)./8))./h;
        advo =  advo - acc.*divv;
        
    elseif strcmp(CTX.SL.Advection(1:3),'FTV')

        velx     =  vel(:,1);
        velz     =  vel(:,2);
        dvx      =  dv(:,1);
        dvz      =  dv(:,2);
        
        u        =  ghost_field(velx,nx,nz,map);
        w        =  ghost_field(velz,nx,nz,map);
        du       =  ghost_field(dvx ,nx,nz,map);
        dw       =  ghost_field(dvz ,nx,nz,map);
        
        velxo    =  oldvel(:,1);
        velzo    =  oldvel(:,2);
        dvxo     =  olddv(:,1);
        dvzo     =  olddv(:,2);
        
        uo       =  ghost_field(velxo,nx,nz,map);
        wo       =  ghost_field(velzo,nx,nz,map);
        duo      =  ghost_field(dvxo ,nx,nz,map);
        dwo      =  ghost_field(dvzo ,nx,nz,map);
        
        wp = (w(jp,ic)+w(jc,ic))./2; 
        wm = (w(jm,ic)+w(jc,ic))./2;
        up = (u(jc,ip)+u(jc,ic))./2;  
        um = (u(jc,im)+u(jc,ic))./2;
        divv = (wp-wm)./h + (up-um)./h;
        
        wp = (w(jp,ic)+dw(jp,ic)+w(jc,ic)+dw(jc,ic))./2;
        wm = (w(jm,ic)+dw(jm,ic)+w(jc,ic)+dw(jc,ic))./2;
        up = (u(jc,ip)+du(jc,ip)+u(jc,ic)+du(jc,ic))./2;  
        um = (u(jc,im)+du(jc,im)+u(jc,ic)+du(jc,ic))./2;
        
        acc = a(jc,ic);
        ajp = a(jp,ic);
        ajm = a(jm,ic);
        aip = a(jc,ip);
        aim = a(jc,im);
        
        adv  =  ((acc+aip)./2.*up - (acc+aim)./2.*um)./h ...
              + ((acc+ajp)./2.*wp - (acc+ajm)./2.*wm)./h;
        adv  =  adv - acc.*divv;
             
        wp = (wo(jp,ic)+wo(jc,ic))./2; 
        wm = (wo(jm,ic)+wo(jc,ic))./2;
        up = (uo(jc,ip)+uo(jc,ic))./2;  
        um = (uo(jc,im)+uo(jc,ic))./2;
        divv = (wp-wm)./h + (up-um)./h;
        
        wp = (wo(jp,ic)+dwo(jp,ic)+wo(jc,ic)+dwo(jc,ic))./2;
        wm = (wo(jm,ic)+dwo(jm,ic)+wo(jc,ic)+dwo(jc,ic))./2;
        up = (uo(jc,ip)+duo(jc,ip)+uo(jc,ic)+duo(jc,ic))./2;  
        um = (uo(jc,im)+duo(jc,im)+uo(jc,ic)+duo(jc,ic))./2;
        
        acc = ao(jc,ic);
        ajp = ao(jp,ic);
        ajm = ao(jm,ic);
        aip = ao(jc,ip);
        aim = ao(jc,im);
        
        advo =  ((acc+aip)./2.*up - (acc+aim)./2.*um)./h ...
              + ((acc+ajp)./2.*wp - (acc+ajm)./2.*wm)./h;
        advo =  advo - acc.*divv;

    else
        adv  = 0;
        advo = 0;
    end

    fieldk(map) =  fieldko(map) - ((adv+advo)./2 - kappa.*(lapl+laplo)./2).*dt;
    fieldk      =  fieldk       + (src+srco)./2.*dt;
    field(:,k)  =  fieldk;
    
end
   
end


function  [fgh]  =  ghost_field(f,nx,nz,map)

ic  =  3:nx+2;
jc  =  3:nz+2;

fgh           =  zeros(nz+4,nx+4);
fgh(jc,ic)    =  f(map);
fgh(   1,ic)  =  f(map(1  ,:));
fgh(   2,ic)  =  f(map(1  ,:));
fgh(nz+3,ic)  =  f(map(end,:));
fgh(nz+4,ic)  =  f(map(end,:));
fgh(jc,   1)  =  f(map(:,  1));
fgh(jc,   2)  =  f(map(:,  1));
fgh(jc,nx+3)  =  f(map(:,end));
fgh(jc,nx+4)  =  f(map(:,end));

end





