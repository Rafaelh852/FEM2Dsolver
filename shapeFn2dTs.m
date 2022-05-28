% this function gives the linear and quadratic shape functions on a
% standard triangle element with vertices (0,0), (0,1) and (1,0)
% [z] = shapeFn2dTs(i,x,y,p)
% input:
%   i indicates the ith shape function psi_i
%   x,y are the input variables
%   p is the order of shape function
% output:
%   z
%
% Last updated: Chung-min Lee April 13, 2017

function [z] = shapeFn2dTs(i,x,y,p)


ind = find( x>= -1e-12 &  y >= -1e-12 & y +x -1 <= 1e-12);
z = zeros(size(x));

if (p == 1)
    switch i
        case 1
            z(ind) = 1 - x(ind) - y(ind);
        case 2
            z(ind) = x(ind);
        case 3
            z(ind) = y(ind);
        otherwise
            display(' wrong input on i')
    end
else % p==2
    switch i
        case 1
            z(ind) = (1 - x(ind) - y(ind)).*(1 - 2*x(ind) - 2*y(ind));
        case 2
            z(ind) = x(ind).*(2*x(ind) - 1);
        case 3
            z(ind) = y(ind).*(2*y(ind) - 1);
        case 4
            z(ind) = 4*x(ind).*(1-x(ind)-y(ind));
        case 5
            z(ind) = 4*x(ind).*y(ind);
        case 6
            z(ind) = 4*y(ind).*(1-x(ind)-y(ind));
        otherwise
            display(' wrong input on i')
    end
end