% this function gives the linear and quadratic shape functions on a triangle
% element with vertices (x1,y1), (x2,y2) and (x3,y3)
% 
% [z] = shapeFn2d(i,x,y,x1,y1,x2,y2,x3,y3,p)
% input:
%   i indicates the ith shape function psi_i
%   x,y are the input variables
%   x1, y1, x2, y2, x3, y3 are coordinates of the vertices
%   p is the order of shape functions
% output:
%   z
%
% Last update: Chung-min Lee April 13, 2017

function [z] = shapeFn2d(i,x,y,x1,y1,x2,y2,x3,y3,p)

XI = ((y3-y1).*(x-x1)+(x1-x3).*(y-y1))./((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
ETA =((y1-y2).*(x-x1)+(x2-x1).*(y-y1))./((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));


[z] = shapeFn2dTs(i,XI,ETA,p);