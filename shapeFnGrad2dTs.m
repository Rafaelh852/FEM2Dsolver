% this program gives the derivatives of the linear and quadratic shape functions on the
% standard triangle element with verivces (0,0), (1,0) and (0,1)
% [psi_xi,psi_eta] = shapeFnGrad2dTs(i,xi,eta,p)
%
% input:
%   i indicates the ith shape function psi_i
%   xi,eta are the input variables
%   p is the order of shape function
% output:
%   psi_xi, psi_eta
%
% Last update: Chung-min Lee April 13, 2017

function [psi_xi,psi_eta] = shapeFnGrad2dTs(i,xi,eta,p)

ind = find( xi>= -1e-12 & eta >= -1e-12 & eta + xi <= 1+1e-12 );
[m,n]=size(xi);
xi = reshape(xi,m*n,1);
eta = reshape(eta,m*n,1);

psi_xi = zeros(m*n,1);
psi_eta = zeros(m*n,1);

if (p == 1)
    switch i
        case 1
            psi_xi(ind) = -1;
            psi_eta(ind) = -1;
        case 2
            psi_xi(ind) = 1;
            psi_eta(ind) = 0;
        case 3
            psi_xi(ind) = 0;
            psi_eta(ind) = 1;
        otherwise
            display(' wrong input on i')
    end
else % p== 2
    
    switch i
        case 1
            psi_xi(ind) = -3 + 4*xi(ind) + 4*eta(ind);
            psi_eta(ind) = -3 + 4*xi(ind) + 4*eta(ind);
        case 2
            psi_xi(ind) = 4*xi(ind) - 1;
            psi_eta(ind) = 0;
        case 3
            psi_xi(ind) = 0;
            psi_eta(ind) = 4*eta(ind) - 1;
        case 4
            psi_xi(ind) = 4 - 8*xi(ind) - 4*eta(ind);
            psi_eta(ind) = -4*xi(ind);
        case 5
            psi_xi(ind) = 4*eta(ind);
            psi_eta(ind) = 4*xi(ind);
        case 6
            psi_xi(ind) = -4*eta(ind);
            psi_eta(ind) = 4 - 4*xi(ind) - 8*eta(ind);
        otherwise
            display(' wrong input on i')
    end
end