%Given: 
    nodes = [-1 0; -0.5 0; -0.5 0.5; 0 -0.5; 0 0; 0 0.5; 0 1; 0.5 0; 0.5 0.5; 1 0];
    triangles = [1 2 3; 2 3 6; 2 5 6; 2 4 5; 4 5 8; 5 6 8; 6 8 9; 8 9 10; 3 6 7; 6 7 9];
    edges = [1 2; 1 3; 2 3; 2 4; 2 5; 2 6; 3 6; 3 7; 4 5; 4 8; 5 6; 5 8; 6 7; 6 8; 6 9; 7 9; 8 9; 8 10; 9 10];
    bdyFn = @(x) 2*(x^2)-0.5;
    bdyFnder = @(x) 4*x;
    k = @(x,y) (x^2)*y+1;

    x = nodes(:,1);
    y = nodes(:,2);

%Project case:
    %f = @(x,y) -4*(5*(x^4)*y + 3*(x^2) + 2*(x^2)*y + 1);
    %y = @(x) 2*(x^2)-0.5;
    %g = @(x,y) 9*(x^4) -4*(x^2) + 11/2;
    % if (y == 0 && abs(x)<=1 && abs(x) >=0.5)
    %     h = @(x,y) 0;
    % elseif (x+y == 1 && x>1 && x<0)
    %     h = @(x,y) (2*sqrt(2))*((x^2)*y+1)*((x^3)+y);
    % elseif (y-x == 1 && x <=0 && x>-1)
    %     h = @(x,y) (2*sqrt(2))*((x^2)*y+1)*((-x^3)+y);
    % end
    
    
%Check case:
    f = @(x,y) -2*x.*y.*pi.*cos(pi*x).*cos(2*pi*y) + (pi^2)*(x.^2*y+1).*sin(pi*x).*cos(2*pi*y) + 2*pi*x.^2.*sin(2*pi*y).*sin(pi*x)+4*pi*pi*(x.^2*y+1).*cos(2*pi*y).*sin(pi*x);
    shapeFn = 1; 
    g = @(x,y) cos(2*pi*y).*sin(pi*x);
    h = @(x,y) fn_h(x,y);
    
%Inputs:
    curveEdge = [0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]';
    bdyNode = [1 1 1 1 0 0 1 1 1 1]';
    bdyEdge = [1 1 0 1 0 0 0 1 0 1 0 0 0 0 0 1 0 1 1]';
    triangleMidPts = [11 12 13; 13 14 15; 14 16 17; 17 18 19; 19 20 21; 16 21 22; 22 23 24; 23 25 26; 15 27 28; 24 28 29];
    midNodes = [-0.25 0.75; -0.75 0.25; -0.25 0.5; -0.5 0.25; -0.75 0; -0.25 0; 0 0.25; 0.25 0.5; 0.25 0.75; 0.5 0.25; 0.75 0.25; 0.75 0; 0.25 -0.25; 0 -0.25; -0.25 -0.25; 0.25 0; 0 0.75; -0.25 0.25];
    Gamma1Nodes = [0 0 0 1 0 0 0 0 0 0]';
    Gamma1MidNodes = zeros(length(midNodes),1);
    Gamma2Edges = bdyEdge - curveEdge;
    Gamma2Nodes = bdyNode - Gamma1Nodes;
    indVec = indVect(Gamma1Nodes,Gamma1MidNodes, shapeFn);

%Refinement: 
    %[nodes,triangles,edges,bdyNode,bdyEdge,curveEdge]=refine(nodes,triangles,edges,bdyNode,bdyEdge,curveEdge,bdyFn,bdyFnder);


%Solving for W:
    [K] = stiffK2d(k, triangles, nodes, triangleMidPts, indVec, shapeFn, 7);
    [F] = loadF2d(k,f,g,h,triangles,nodes,edges,triangleMidPts,midNodes,indVec,Gamma1Nodes,Gamma1MidNodes,Gamma2Edges,shapeFn,7,3);
    [Kbar,Q] = permRCM(K); %Kbar = Q*K*(Q)'
    Fbar = Q*F;

    p = bandwidth(Kbar,'lower');
    q = bandwidth(Kbar,'upper');
    LU = bandLU(Kbar,p,q);
    L = tril(LU,-1) + diag(ones(1,length(LU)));
    U = triu(LU);

    y = bandforward(L,Fbar,p);
    Wbar = bandbackward(U,y,p);
    

    W = Q.'*Wbar;

function [val] = fn_h(x,y)
    val = zeros(length(x),1);
        
    ind1 = ((abs(y) < 1e-10) & (abs(x) > 0.5-1e-10) & (abs(x) < 1+1e-10));
    val(ind1) = 0;
        
    ind2 = ((abs(x+y-1) < 1e-10) & (x > 0) & (x < 1));
    val(ind2) = pi/sqrt(2)*(cos(pi*x(ind2)).*cos(2*pi*y(ind2))-2*sin(pi*x(ind2)).*sin(2*pi*y(ind2)));
        
    ind3 = ((abs(y-x-1) < 1e-10) & (x > -1) & (x < 1e-10));
    val(ind3) = -pi/sqrt(2)*(cos(pi*x(ind3)).*cos(2*pi*y(ind3))+2*sin(pi*x(ind3)).*sin(2*pi*y(ind3)));
end
