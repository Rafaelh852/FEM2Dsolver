% this function takes the computed values of the nodes and gives the
% approximated solution on the domain
% 
% [uh] =
% approxSol2d(w,g,indVec,nodes,triangles,midNodes,triangleMidPts,shapeFn)
% 
% where w is the values of the nodes, 
% g is the Dirichlet BC on Gamma_1, 
% indVec is the index array maps original node numbering to free node
% numbering.  
% nodes are the corrdinates of nodes, 
% triangles records the vertices of the triangles, 
% triangleMidPts records the node number of midedge points, 
% midNodes has the coordinates of the midedge points, 
% shapeFn indicates the choice of linear or quadratic shape functions 
% Last updated: Chung-min Lee April 24, 2021

function [uh] = approxSol2d(w,g,indVec,nodes,triangles,midNodes,triangleMidPts,shapeFn)

if (shapeFn == 1) % if linear shape functions are used
    uh = @(x,y) linearApprox2d(x,y,w,g,indVec,nodes,triangles);
else    % if quadratic shape functions are used
    uh = @(x,y) quadraticApprox2d(x,y,w,g,indVec,nodes,triangles,midNodes,triangleMidPts);
end


%%---------
% assemble linear shape functions on the whole domain

function z = linearApprox2d(x,y,w,g,indVec,nodes,triangles)

noNodes = length(nodes(:,1));

% put all node values together
uu = zeros(noNodes,1);
uu(logical(indVec)) = w;
gamma1s = (indVec==0);
uu(gamma1s) = g(nodes(gamma1s,1),nodes(gamma1s,2));

N1 = triangles(:,1);
N2 = triangles(:,2);
N3 = triangles(:,3);

notri = length(N1);

[lxr,lxc] = size(x);
lx = lxr*lxc;
z = zeros(lx,1);
if (lx > notri)
    ckNan = zeros(lx,1);
    for i = 1:notri
        psi1 = shapeFn2d(1,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),1);
        psi2 = shapeFn2d(2,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),1);
        psi3 = shapeFn2d(3,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),1);
        ckNan = ckNan+(psi1 | psi2 | psi3);
        z = z + uu(N1(i)).*psi1 + uu(N2(i)).*psi2 + uu(N3(i)).*psi3;
    end
    z(logical(ckNan)) = z(logical(ckNan))./ckNan(logical(ckNan));
    z(~ckNan) = nan;
else
     for i = 1:lx
         psi1 = shapeFn2d(1,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),1);
         psi2 = shapeFn2d(2,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),1);
         psi3 = shapeFn2d(3,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),1);
         ckNan = sum((psi1 | psi2 | psi3),1);
         if (ckNan)
            z(i) = sum(uu(N1).*psi1 + uu(N2).*psi2 + uu(N3).*psi3,1)/ckNan;
         else
             z(i) = nan;
         end
     end
end
z = reshape(z,lxr,lxc);

%%------------------
% assemble quadratic shape functions on the whole domain

function z = quadraticApprox2d(x,y,w,g,indVec,nodes,triangles,midNodes,triangleMidPts)

noVer = length(nodes(:,1));
noNodes = noVer +length(midNodes(:,1));

% put all node values together
uu = zeros(noNodes,1);
uu(logical(indVec)) = w;
gamma1s = find(indVec==0);
gamma1nodes = find(gamma1s <= noVer);
gamma1midNodes = setdiff([1:length(gamma1s)],gamma1nodes);
uu(gamma1s(gamma1nodes)) = g(nodes(gamma1s(gamma1nodes),1),nodes(gamma1s(gamma1nodes),2));
uu(gamma1s(gamma1midNodes)) = g(midNodes(gamma1s(gamma1midNodes)-noVer,1),midNodes(gamma1s(gamma1midNodes)-noVer,2));

%all vertices and mid-edge nodes
N1 = triangles(:,1);
N2 = triangles(:,2);
N3 = triangles(:,3);
N4 = triangleMidPts(:,1);
N5 = triangleMidPts(:,2);
N6 = triangleMidPts(:,3);

notri = length(N1);

[lxr,lxc] = size(x);
lx = lxr*lxc;
z = zeros(lx,1);
if (lx > notri)
    ckNan = zeros(lx,1);
    for i = 1:notri
        psi1 = shapeFn2d(1,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2);
        psi2 = shapeFn2d(2,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2);
        psi3 = shapeFn2d(3,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2);
        psi4 = shapeFn2d(4,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2);
        psi5 = shapeFn2d(5,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2);
        psi6 = shapeFn2d(6,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2);
        ckNan = ckNan+(psi1 | psi2 | psi3 | psi4 | psi5 | psi6);
        z = z+uu(N1(i)).*psi1 + uu(N2(i)).*psi2+ uu(N3(i)).*psi3 + uu(N4(i)).*psi4 + uu(N5(i)).*psi5 + uu(N6(i)).*psi6;
    end
    z(logical(ckNan)) = z(logical(ckNan))./ckNan(logical(ckNan));
    z(~ckNan) = nan;
else
    for i = 1:lx
        psi1 = shapeFn2d(1,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2);
        psi2 = shapeFn2d(2,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2);
        psi3 = shapeFn2d(3,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2);
        psi4 = shapeFn2d(4,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2);
        psi5 = shapeFn2d(5,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2);
        psi6 = shapeFn2d(6,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2);
        ckNan = sum((psi1 | psi2 | psi3 | psi4 | psi5 |psi6),1);
        if (ckNan)
            z(i) = sum(uu(N1).*psi1 + uu(N2).*psi2 + uu(N3).*psi3 + uu(N4).*psi4 + uu(N5).*psi5 + uu(N6).*psi6,1)/ckNan;
        else
            z(i) = nan;
        end
    end
end
z = reshape(z,lxr,lxc);


    
