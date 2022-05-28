function [indVec] = indVect(Gamma1Nodes,Gamma1MidNodes, shapeFn)
k = 0; 

if shapeFn==2 % if quadratic, we want to check the midnodes too. 
    Gamma1Nodes = [ Gamma1Nodes;Gamma1MidNodes];
end

indVec = zeros(length(Gamma1Nodes),1);

for i = 1: length(indVec)
   if ( Gamma1Nodes(i) == 1 )% node is on gamma1 
       indVec(i) = 0;% constrained node 
   else
       k = k + 1;
       indVec(i) = k; % free node number
   end
   
end


end