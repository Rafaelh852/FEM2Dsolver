% this program gives the derivatives of the linear and quadratic shape functions on a 
% triangle element with verivces (x1,y1), (x2,y2) and (x3,y3)
% [gradpsi] = shapeFnGrad2d(i,x,y,x1,y1,x2,y2,x3,y3,p)
% 
% input:
%   i indicates the ith shape function psi_i
%   x,y are the input variables
%   x1, y1, x2, y2, x3, y3 are coordinates of vertices 
%   p is the order of shape function
% output:
%   gradpsi
%
% Last update: Chung-min Lee April 16, 2021

function [gradpsi] = shapeFnGrad2d(i,x,y,x1,y1,x2,y2,x3,y3,p)

xi = ((y3-y1).*(x-x1)+(x1-x3).*(y-y1))./((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
eta =((y1-y2).*(x-x1)+(x2-x1).*(y-y1))./((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));

[a,b]=shapeFnGrad2dTs(i,xi,eta,p);

psi_x =( ((y3-y1).*a + (y1-y2).*b)./((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)));
psi_y = (((x1-x3).*a + (x2-x1).*b)./((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)));
[gradpsi] = [psi_x; psi_y];

