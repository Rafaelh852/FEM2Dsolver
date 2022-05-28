function [z] = bandforward(L,f,p)

n = length(f);

for j = 1:n
    for i= j+1:min([j+p,n])
        f(i) = f(i) - L(i,j)*f(j);
    end 
end
 z = f;
end