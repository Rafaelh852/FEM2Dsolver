clear
nodes = [-1 0; -0.5 0; -0.5 0.5; 0 -0.5; 0 0; 0 0.5; 0 1; 0.5 0; 0.5 0.5; 1 0];
triangles = [1 2 3; 2 3 6; 2 5 6; 2 4 5; 4 5 8; 5 6 8; 6 8 9; 8 9 10; 3 6 7; 6 7 9];
edges = [ 1 2 ; 1 3; 2 3; 2 4; 2 5; 2 6; 3 6; 3 7; 4 5; 4 8; 5 6; 5 8; 6 7; 6 8; 6 9; 7 9; 8 9; 8 10; 9 10];


bdyFn = @(x)(2*x.^2 - .5);
bdyFnder = @(x)(4*x);
% computer manually
bdyNode = [1 1 1 1 0 0 1 1 1 1]'; 
bdyEdge = [1;1;0;1;0;0;0;1;0;1;0;0;0;0;0;1;0;1;1]; 
curveEdge = [0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]';

ut = @(x,y)(x.^4+2*y.^2+5);
utg = @(x,y)([4*x.^3; 4*y]);
% proj 2 data
k = @(x,y)(x.^2.*y + 1);
f = @(x,y)(-4*(5*x.^4.*y + 3*x.^2 + 2*x.^2.*y + 1));
g = @(x,y)(9*x.^4 - 4*x.^2 + 11/2);
h = @(x,y)(fn(x,y));

%"check" data bellow

% k = @(x,y) x.^2*y+1;
% f = @(x,y) -2*x.*y.*pi.*cos(pi*x).*cos(2*pi*y) + pi^2*(x.^2*y+1).*sin(pi*x).*cos(2*pi*y)+2*pi*x.^2.*sin(2*pi*y).*sin(pi*x)+4*pi*pi*(x.^2*y+1).*cos(2*pi*y).*sin(pi*x);
% g = @(x,y) cos(2*pi*y).*sin(pi*x); 
% h = @(x,y) fn_h(x,y);

shapeFn = 1;
maxh = 0.8;

x = nodes(:,1);
y = nodes(:,2);

xx = -1:0.02:1;
yy = -0.5:0.02:1;
[XX,YY] = meshgrid(xx,yy);

format shortg;
[uh, graduh, hMax, nodes1, triangles1] = myFE2dbvp(k, f, g, h, maxh, nodes, triangles, edges, bdyNode, bdyEdge, curveEdge, bdyFn, bdyFnder, shapeFn);

x1 = [-0.74 -0.3 0.16 0.8]';

y1 = [0.1 0.54 -0.23 0.09]';

uh(x1,y1)


% H1norm2d(uh,graduh,nodes1,triangles1,ut,utg)
subplot(3,1,1)
surf(XX,YY,uh(XX,YY));
title("abs. error for linear shape functions with maxh =0.2");
xlabel("x");
ylabel("y");
zlabel("|u_{true}-uh|");

shapeFn = 2;
maxh = 0.4;
[uh, graduh, hMax, nodes2, triangles2] = myFE2dbvp(k, f, g, h, maxh, nodes, triangles, edges, bdyNode, bdyEdge, curveEdge, bdyFn, bdyFnder, shapeFn);

subplot(4,1,2)
surf(XX,YY,uh(XX,YY));
title("abs. error for quadratic shape functions with maxh =0.2");
xlabel("x");
ylabel("y");
zlabel("|u_{true}-uh|");

 maxh = 0.2;
[uh, graduh, hMax, nodes3, triangles3] = myFE2dbvp(k, f, g, h, maxh, nodes, triangles, edges, bdyNode, bdyEdge, curveEdge, bdyFn, bdyFnder, shapeFn);

subplot(4,1,3)
surf(XX,YY,(uh(XX,YY)));
title("abs. error for maxh =0.2");
xlabel("x");
ylabel("y");
zlabel("|u_{true}-uh|");

subplot(4,1,4)
surf(XX,YY,(ut(XX,YY)));
title("abs. error for maxh =0.2");
xlabel("x");
ylabel("y");
zlabel("|u_{true}-uh|");


function [val] = fn(x,y)
    val = zeros(length(x),1);
        
    ind1 = ((abs(y) < 1e-10) & (abs(x) > 0.5-1e-10) & (abs(x) < 1+1e-10));
    val(ind1) = 0;
        
    ind2 = ((abs(x+y-1) < 1e-10) & (x > 0) & (x < 1));
    val(ind2) = 2*sqrt(2)*(x(ind2).^2.*y(ind2) + 1).*(x(ind2).^3 + y(ind2))  ;
        
    ind3 = ((abs(y-x-1) < 1e-10) & (x > -1) & (x < 1e-10));
    val(ind3) =2*sqrt(2)*(x(ind3).^2.*y(ind3) + 1).*(-x(ind3).^3+ y(ind3)) ;

end

function [val] = fn_h(x,y)
    val = zeros(length(x),1);
        
    ind1 = ((abs(y) < 1e-10) & (abs(x) > 0.5-1e-10) & (abs(x) < 1+1e-10));
    val(ind1) = 0;
        
    ind2 = ((abs(x+y-1) < 1e-10) & (x > 0) & (x < 1));
    val(ind2) = pi/sqrt(2)*(cos(pi*x(ind2)).*cos(2*pi*y(ind2))-2*sin(pi*x(ind2)).*sin(2*pi*y(ind2)));
        
    ind3 = ((abs(y-x-1) < 1e-10) & (x > -1) & (x < 1e-10));
    val(ind3) = -pi/sqrt(2)*(cos(pi*x(ind3)).*cos(2*pi*y(ind3))+2*sin(pi*x(ind3)).*sin(2*pi*y(ind3)));
end