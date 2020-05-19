% calculate value of shape function and shape function derivatives at given
% coordinates
%
% [N,dNdS] = ShapeFuncts(S,type)
% Input is an array S of dimensions (2,n), containing coordinates of points
% (typically integration points or lagrangian particles), at which shape 
% functions and derivatives will be calculated.
% Function returns shape functions in array N (nodes_el,n), and shape
% function derivatives in array dNdS (2,nodes_el,n)
% linear or quadratic elements supported
%
% sorting of nodes in element is
% 
% type = P0:    .---.  
%               | o | 1
%               .---.  
%
%
% type = Q1:  1 o---o 3
%               |   |
%             2 o---o 4
%
%                   4
% type = Q2:  1 o---o---o 7
%               |   |5  |
%             2 o---o---o 8
%               |   |   |
%             3 o---o---o 9
%                   6   

% created  20140729 Tobias Keller
% modified 20190418 Tobias Keller


function  [N,dNdS]  =  ShapeFuncts(S,type)


%***  P0 - piece-wise constant shape functions

if strcmp(type,'P0')
    
    n             =  length(S(1,:));
    N             =  zeros(1,n);
    dNdS          =  zeros(2,1,n);
 
    z             =  S(1,:);
    x             =  S(2,:);
    
    N(1,:)        =  1;
    
    dNdS(1,1,:)   =  1;
    dNdS(2,1,:)   =  1;
    
    
%***  Q1 - linear shape functions

elseif strcmp(type,'Q1')
    
    n             =  length(S(1,:));
    N             =  zeros(4,n);
    dNdS          =  zeros(2,4,n);
 
    z             =  S(1,:);
    x             =  S(2,:);
    
    N(1,:)        =  0.25 .* (1-x) .* (1-z);
    N(2,:)        =  0.25 .* (1-x) .* (1+z);
    N(3,:)        =  0.25 .* (1+x) .* (1-z);
    N(4,:)        =  0.25 .* (1+x) .* (1+z);
    
    dNdS(1,1,:)   =  0.25.*(-1+z);
    dNdS(1,2,:)   =  0.25.*(-1-z);
    dNdS(1,3,:)   =  0.25.*(+1-z);
    dNdS(1,4,:)   =  0.25.*(+1+z);
    
    dNdS(2,1,:)   =  0.25.*(-1+x);
    dNdS(2,2,:)   =  0.25.*(+1-x);
    dNdS(2,3,:)   =  0.25.*(-1-x);
    dNdS(2,4,:)   =  0.25.*(+1+x);
    
    
%***  Q2 - quadratic shape functions

elseif strcmp(type,'Q2')
    
    n            =  length(S(1,:));
    N            =  zeros(9,n);
    dNdS         =  zeros(2,9,n);
    
    z            =  S(1,:);
    x            =  S(2,:);

    N(1,:)       =  (x.^2 - x) .* (z.^2 - z) .* 0.25;
    N(2,:)       =  (x.^2 - x) .* (1 - z.^2) .* 0.50;
    N(3,:)       =  (x.^2 - x) .* (z.^2 + z) .* 0.25;

    N(4,:)       =  (1 - x.^2) .* (z.^2 - z) .* 0.50;
    N(5,:)       =  (1 - x.^2) .* (1 - z.^2) .* 1.00;
    N(6,:)       =  (1 - x.^2) .* (z.^2 + z) .* 0.50;

    N(7,:)       =  (x.^2 + x) .* (z.^2 - z) .* 0.25;
    N(8,:)       =  (x.^2 + x) .* (1 - z.^2) .* 0.50;
    N(9,:)       =  (x.^2 + x) .* (z.^2 + z) .* 0.25;
    
    
    dNdS(1,1,:)  =  (2.*x - 1) .* (z.^2 - z) .* 0.25;
    dNdS(1,2,:)  =  (2.*x - 1) .* (1 - z.^2) .* 0.50;
    dNdS(1,3,:)  =  (2.*x - 1) .* (z.^2 + z) .* 0.25;

    dNdS(1,4,:)  =  (  - 2.*x) .* (z.^2 - z) .* 0.50;
    dNdS(1,5,:)  =  (  - 2.*x) .* (1 - z.^2) .* 1.00;
    dNdS(1,6,:)  =  (  - 2.*x) .* (z.^2 + z) .* 0.50;

    dNdS(1,7,:)  =  (2.*x + 1) .* (z.^2 - z) .* 0.25;
    dNdS(1,8,:)  =  (2.*x + 1) .* (1 - z.^2) .* 0.50;
    dNdS(1,9,:)  =  (2.*x + 1) .* (z.^2 + z) .* 0.25;

    dNdS(2,1,:)  =  (x.^2 - x) .* (2.*z - 1) .* 0.25;
    dNdS(2,2,:)  =  (x.^2 - x) .* (  - 2.*z) .* 0.50;
    dNdS(2,3,:)  =  (x.^2 - x) .* (2.*z + 1) .* 0.25;

    dNdS(2,4,:)  =  (1 - x.^2) .* (2.*z - 1) .* 0.50;
    dNdS(2,5,:)  =  (1 - x.^2) .* (  - 2.*z) .* 1.00;
    dNdS(2,6,:)  =  (1 - x.^2) .* (2.*z + 1) .* 0.50;

    dNdS(2,7,:)  =  (x.^2 + x) .* (2.*z - 1) .* 0.25;
    dNdS(2,8,:)  =  (x.^2 + x) .* (  - 2.*z) .* 0.50;
    dNdS(2,9,:)  =  (x.^2 + x) .* (2.*z + 1) .* 0.25;
    
end


