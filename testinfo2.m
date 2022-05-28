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
gradut = @(x,y)([4*x.^3; 4*y]);
% proj 2 data
k = @(x,y)(x.^2.*y + 1);
f = @(x,y)(-4*(5*x.^4.*y + 3*x.^2 + 2*x.^2.*y + 1));
g = @(x,y)(9*x.^4 - 4*x.^2 + 11/2);
h = @(x,y)(fn(x,y));

format shortg;

shapeFn = 1;
maxh = 0.8;

[uh, graduh, hMax, nodes1, triangles1] = myFE2dbvp(k, f, g, h, maxh, nodes, triangles, edges, bdyNode, bdyEdge, curveEdge, bdyFn, bdyFnder, shapeFn);

udiff = @(x,y)(abs(ut(x,y)-uh(x,y)));
%graddiff = @(x,y)(gradut(x,y)-graduh(x,y)');

%[L2norm,H1norm] = H1norm2d(udiff,graddiff,udiff,nodes1,triangles1)


L2normVector1 = [0.38759; 0.082026; 0.026851];
L2normVector2 = [0.14074; 0.021789; 0.005427];

H1normVector1 = [0.94353; 0.42248; 0.21168];
H1normVector2 = [0.21151; 0.13502; 0.055603];

H1 = [0.70711; 0.369; 0.18482];
H2 = [ 0.70711; 0.369; 0.18482];

logH = [ones(3,1), log(H1)];
lognorm = log(H1normVector2);
rate_p = logH\lognorm;
rate_p(2)



function [val] = fn(x,y)
    val = zeros(length(x),1);
        
    ind1 = ((abs(y) < 1e-10) & (abs(x) > 0.5-1e-10) & (abs(x) < 1+1e-10));
    val(ind1) = 0;
        
    ind2 = ((abs(x+y-1) < 1e-10) & (x > 0) & (x < 1));
    val(ind2) = 2*sqrt(2)*(x(ind2).^2.*y(ind2) + 1).*(x(ind2).^3 + y(ind2))  ;
        
    ind3 = ((abs(y-x-1) < 1e-10) & (x > -1) & (x < 1e-10));
    val(ind3) =2*sqrt(2)*(x(ind3).^2.*y(ind3) + 1).*(-x(ind3).^3+ y(ind3)) ;

end