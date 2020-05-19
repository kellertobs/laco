
function  [] = PlotField(field,FE,style)

load FireAndIce.mat
if FE.nxEl == 1; style = 'one'; end

k = size(field,2);

if     length(field) == FE.NEl
    type = 'El';
elseif length(field) == FE.NQ1
    type = 'Q1';
elseif length(field) == FE.NQ2
    type = 'Q2';
elseif length(field) == FE.NIP
    type = 'IP';
else
    type = 'unknown';
end

if FE.D > 1e3
    scale = 1e3;
else
    scale = 1;
end

switch type
    
    case 'El'
        X = reshape(FE.CoordEl(:,1)./scale,FE.nzEl,FE.nxEl);
        Z = reshape(FE.CoordEl(:,2)./scale,FE.nzEl,FE.nxEl);
        A = reshape(field                 ,FE.nzEl,FE.nxEl,k);
    
    case 'Q1'
        X = reshape(FE.CoordQ1(:,1)./scale,FE.nzQ1,FE.nxQ1);
        Z = reshape(FE.CoordQ1(:,2)./scale,FE.nzQ1,FE.nxQ1);
        A = reshape(field                 ,FE.nzQ1,FE.nxQ1,k);
        
    case 'Q2'
        X = reshape(FE.CoordQ2(:,1)./scale,FE.nzQ2,FE.nxQ2);
        Z = reshape(FE.CoordQ2(:,2)./scale,FE.nzQ2,FE.nxQ2);
        A = reshape(field                 ,FE.nzQ2,FE.nxQ2,k);
        
    case 'IP'
        X = reshape(FE.CoordIP(:,1)./scale,FE.nzIP,FE.nxIP);
        Z = reshape(FE.CoordIP(:,2)./scale,FE.nzIP,FE.nxIP);
        A = reshape(field                 ,FE.nzIP,FE.nxIP,k);
        
    otherwise
        error('The input data has an unknown format! Check and try again.');
end

switch style(end-2:end)
    
    case '_rr'
        X = [X,2.*FE.W./scale-fliplr(X)];
        Z = [Z,fliplr(Z)];
        A = [A,fliplr(A)];

    case '_rl'
        X = [FE.W./scale-fliplr(X),X+FE.W./scale];
        Z = [fliplr(Z),Z];
        A = [fliplr(A),A];
        
    case 'rrU'
        X = [X,2.*FE.W./scale-fliplr(X)];
        Z = [Z, fliplr(Z)];
        A = [A,-fliplr(A)];

    case 'rlU'
        X = [FE.W./scale-fliplr(X),X+FE.W./scale];
        Z = [ fliplr(Z),Z];
        A = [-fliplr(A),A];
end

switch style(1:3)
    
    case 'img'
        imagesc(mean(X,1),mean(Z,2),A);
        shading flat;
        axis ij equal tight; 
        colormap(FireAndIce);
        colorbar('TicklabelInterpreter','Latex','FontSize',15);

    case 'srf'
        surf(X,Z,A);
        view(0,90);
        shading flat;
        axis ij equal tight; 
        colormap(FireAndIce);
        colorbar('TicklabelInterpreter','Latex','FontSize',15);
        
     case 'qui'
        imagesc(mean(X,1),mean(Z,2),A(:,:,3)); hold on;
        quiver(X(2:6:end-1,2:6:end-1),Z(2:6:end-1,2:6:end-1),A(2:6:end-1,2:6:end-1,1),A(2:6:end-1,2:6:end-1,2),2,'w','LineWidth',1.5);
        shading flat;
        axis ij equal tight; 
        colormap(FireAndIce);
        colorbar('TicklabelInterpreter','Latex','FontSize',15);
        
    case 'one'
        plot(mean(A,2),mean(Z,2));
        axis ij tight;
        
    otherwise
        error('Enter a valid plot style: img (image), srf (surface), qui (quiver), one (1-D line plot)');
        
end

box on;

end