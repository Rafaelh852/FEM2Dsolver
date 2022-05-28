function hmax = findhmax(edges, nodes)
% here we make a vector full of the lengths of the edges. 
edgelengths = zeros( length(edges(:,1)),1);

for i = 1: length(edges(:,1))
    
    x1 = nodes(edges(i,1),1); 
    y1 = nodes(edges(i,1),2);
    
    x2 = nodes(edges(i,2),1);
    y2 = nodes(edges(i,2),2);
     
    edgelengths(i) = sqrt( (x2 - x1)^2 + (y2-y1)^2);    
end

hmax = max(edgelengths);

end