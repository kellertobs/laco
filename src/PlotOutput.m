% preview MOR model output as time series

function  [CTX] = PlotOutput(CTX,start,step,stop,name,varargin)

load FireAndIce;

holdfig  =  0;
flipcmap =  0;
scale    = 'lin';
oper     = 'none';
color    = 'r';
comp     =  1;
style    =  'img';
hzv      =  false;
printfig =  false;

format  =  '-dpng';
res     =  '-r300';
rend    =  '-opengl';
    
n = 0;
while n < length(varargin)
    n = n+1;
    if ischar(varargin{n})
        if (strcmp(varargin{n}(1:3),'sub') || strcmp(varargin{n}(1:3),'add') ...
         || strcmp(varargin{n}(1:3),'div') || strcmp(varargin{n}(1:3),'mul'))
            oper = varargin{n};
            n    = n+1;
            name2 = varargin{n};
        elseif strcmp(varargin{n}(1:3),'fli')
            flipcmap = 1;
        elseif strcmp(varargin{n}(1:3),'hol')
            holdfig = 1;
        elseif (strcmp(varargin{n}(1:3),'col') || strcmp(varargin{n}(1:3),'col'))
            n = n+1;
            color = varargin{n};
        elseif (strcmp(varargin{n}(1:3),'com'))
            n = n+1;
            comp = varargin{n};
        elseif (strcmp(varargin{n}(1:3),'hzv'))
            hzv = true;
        elseif (strcmp(varargin{n}(1:3),'pri'))
            printfig = true;
        elseif (strcmp(varargin{n}(1:3),'lin') || strcmp(varargin{n}(1:3),'log') || strcmp(varargin{n}(1:3),'abs'))
            scale = varargin{n};
        else
            style = varargin{n};
        end
    elseif isvector(varargin{n})
        clim = varargin{n};
    end
end

if strcmp(style(1:3),'qui') || strcmp(name(1),'V')
    if (strcmp(name,'U') || strcmp(name,'W') || strcmp(name,'V'))
        name  = 'V';
        uname = 'U';
        wname = 'W';
    elseif length(name)==2 && (strcmp(name(1),'U') || strcmp(name(1),'W') || strcmp(name(1),'V'))
        name  = ['V',name(2)];
        uname = ['U',name(2)];
        wname = ['W',name(2)];
    else
        disp('Need to pass a velocity field for quiver plot!');
        disp('Ignoring quiver, plotting as image instead.');
        scale = 'lin';
    end
end

fh = figure(100);
if holdfig
    hold on;
else
    clf;
end
    
for frame=start:step:stop
    
    if strcmp(style(1:3),'qui') || strcmp(name(1),'V')
        CTX.SL.(name) = sqrt(CTX.SL.(uname).^2 + CTX.SL.(wname).^2);
    end
    
    time = CTX.TIME.total;
    mn = 60;
    hr = 3600;
    dy = 3600*24;
    yr = 3600*24*365.25;
    if time < mn
        tunit = ' sec';
    elseif time >= mn && time < hr
        tunit = ' min';
        time  = time/mn;
    elseif time >=  hr && time < dy
        tunit = ' hr';
        time  = time/hr;
    elseif time >=  dy && time < yr
        tunit = ' day';
        time  = time/dy;
    elseif time >=  yr && time < 1e3*yr
        tunit = ' yr';
        time  = time/yr;
    elseif time >=  1000*yr && time < 1e6*yr
        tunit = ' kyr';
        time  = time/1e3/yr;
    else
        tunit = ' Myr';
        time  = time/1e6/yr;
    end
    
    if any(strcmp(fieldnames(CTX.SL),name))
        A = CTX.SL;
    elseif any(strcmp(fieldnames(CTX.MP),name))
        A = CTX.MP;
    else
        error('Specified field not found in any of the output structures.')
    end

    if strcmp(style(1:3),'qui')
        field      = zeros(length(CTX.SL.(uname)),3);
        field(:,1) = A.(uname);
        field(:,2) = A.(wname);
        field(:,3) = A.( name);
    else
        field  = A.(name)(:,comp);
        if exist('name2','var') && ischar(name2); field2 = A.(name2)(:,comp); end
    end
    
    if strcmp(oper(1:3),'sub')
        if ischar(name2)
            field = field - field2;
        else
            field = field - name2;
        end
    elseif strcmp(oper(1:3),'add')
        if ischar(name2)
            field = field + field2;
        else
            field = field + name2;
        end
    elseif strcmp(oper(1:3),'div')
        if ischar(name2)
            field = field ./ field2;
        else
            field = field ./ name2;
        end
    elseif strcmp(oper(1:3),'mul')
        if ischar(name2)
            field = field .* field2;
        else
            field = field .* name2;
        end
    end
        
    if strcmp(scale(1:3),'log'); field = log10(abs(field)); end
    if strcmp(scale(1:3),'abs'); field = abs(field); end
    
    if hzv
        if any(strcmp(name,{'U';'W';'V'}))
            k      =  size(field,2);
            wk     =  reshape(field,CTX.FE.nzQ2,CTX.FE.nxQ2,k);
            wk     =  wk - mean(wk,2);
            field  =  reshape(wk,CTX.FE.NQ2,k);
        else
            k      =  size(field,2);
            wk     =  reshape(field,CTX.FE.nzQ1,CTX.FE.nxQ1,k);
            wk     =  wk - mean(wk,2);
            field  =  reshape(wk,CTX.FE.NQ1,k);
        end
    end
    
    if CTX.FE.nxEl<=2
        if CTX.FE.nzEl<=2
            style = 'zero';
        else
            style = 'one';
        end
    end            
    
    if CTX.FE.D > 1e3
        xunit = 'km';
    else
        xunit = 'm';
    end

    if strcmp(style,'zero')
        
        plot(time,mean(field),'o','Color',color,'MarkerSize',10,'LineWidth',2);
        
        xlabel(['Time ' tunit],'Interpreter','Latex','FontSize',18);

    else
        
        PlotField(field,CTX.FE,style);
        
        if flipcmap
            colormap(flipud(colormap));
        end
        
        if exist('clim','var')
            if strcmp(style,'one')
                set(gca,'XLim',clim);
                set(gca,'YLim',[0,CTX.FE.D]);
            else
                caxis(clim)
            end
        end
        
        xlabel(['Width [' xunit ']'],'Interpreter','Latex','FontSize',18);
        ylabel(['Depth [' xunit ']'],'Interpreter','Latex','FontSize',18);
        
    end
    
    set(gca,'TicklabelInterpreter','Latex','FontSize',15);
    title([name,'; ',num2str(time,'%1.2f'),tunit],'Interpreter','Latex','FontSize',18); drawnow;
    
    if printfig
        print(fh,format,res,rend,[CTX.IO.DataDir '/' CTX.IO.RunID '/' CTX.IO.RunID '_' name '_' num2str(frame)],'-loose');
    end
end

if holdfig
    hold off;
end

    
end




