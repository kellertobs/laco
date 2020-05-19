% LivePlotting    EDIFICE: Plot model results live during run
%
% []  =  LivePlotting(CTX)
%
%   Function creates figures and plots model solution and material point 
%   properties live during model run. 
%
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function  []  =  LivePlotting(CTX)

if strcmp(CTX.IO.LivePlot,'ON')
    
    if CTX.FE.nxEl >= CTX.FE.nzEl
        sp1 = 211;
        sp2 = 212;
    else
        sp1 = 121;
        sp2 = 122;
    end
    
    n=1;
    
    if std(CTX.MP.Mat)>=1e-6
        figure(n); n=n+1; clf;
        subplot(sp1);
        PlotField(CTX.MP.Mat,CTX.FE,CTX.IO.PlotStyle);
        title('Material types')
        subplot(sp2);
        PlotField(CTX.MP.Rho,CTX.FE,CTX.IO.PlotStyle);
        title('Density [kg/m3]')
        drawnow
    end
    
    figure(n); n=n+1; clf;
    subplot(sp1);
    PlotField( CTX.SL.U*CTX.TIME.spyr*1000,CTX.FE,[CTX.IO.PlotStyle,'U']);
    title('x-Velocity [mm/yr]')
    subplot(sp2);
    PlotField(-CTX.SL.W*CTX.TIME.spyr*1000,CTX.FE,CTX.IO.PlotStyle);
    title('z-Velocity [mm/yr]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(311);
    PlotField(CTX.MP.DivV,CTX.FE,CTX.IO.PlotStyle);
    title('Volume strain rate [log10 1/s]')
    subplot(312);
    PlotField(log10(CTX.MP.EII(:,1)),CTX.FE,CTX.IO.PlotStyle);
    title('Shear strain rate [log10 1/s]')
    subplot(313);
    PlotField(log10(CTX.MP.TII),CTX.FE,CTX.IO.PlotStyle);
    title('Shear stress [log10 Pa]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(sp1);
    PlotField(CTX.SL.P./1e6,CTX.FE,CTX.IO.PlotStyle);
    title('Dynamic pressure [MPa]')
    subplot(sp2);
    PlotField(CTX.SL.Pt./1e6,CTX.FE,CTX.IO.PlotStyle);
    title('Total pressure [MPa]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(2,2,1);
    PlotField(log10(CTX.MP.Eta),CTX.FE,CTX.IO.PlotStyle);
    title('Eta V [Pas]')
    subplot(2,2,2);
    PlotField(log10(CTX.MP.EtaVP),CTX.FE,CTX.IO.PlotStyle);
    title('Eta VP [Pas]')
    subplot(2,2,3);
    PlotField((CTX.MP.Chi),CTX.FE,CTX.IO.PlotStyle);
    title('Chi VEP [1]')
    subplot(2,2,4);
    PlotField(log10(CTX.MP.EtaVEP),CTX.FE,CTX.IO.PlotStyle);
    title('Eta VEP [Pas]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(sp1);
    PlotField(log10(CTX.MP.YieldStr),CTX.FE,CTX.IO.PlotStyle);
    title('Yield Stress [log10 Pa]')
    subplot(sp2);
    PlotField(log10(max(CTX.RHEO.Dmg0.*1e-3,CTX.MP.Dmg)),CTX.FE,CTX.IO.PlotStyle);
    title('Damage strain [log10]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(2,2,1);
    PlotField(max(-18,log10(CTX.MP.EII(:,1))),CTX.FE,CTX.IO.PlotStyle);
    title('Dev. strainr total [1/s]')
    subplot(2,2,2);
    PlotField(max(-18,log10(CTX.MP.EII(:,2))),CTX.FE,CTX.IO.PlotStyle);
    title('Dev. strainr viscous [1/s]')
    subplot(2,2,3);
    PlotField(max(-18,log10(CTX.MP.EII(:,3))),CTX.FE,CTX.IO.PlotStyle);
    title('Dev. strainr elastic [1/s]')
    subplot(2,2,4);
    PlotField(max(-18,log10(CTX.MP.EII(:,4))),CTX.FE,CTX.IO.PlotStyle);
    title('Dev. strainr plastic [1/s]')
    drawnow
    
end

end