function [z] = bandbackward(L,f,p)

n = length(f);
  
for j = n:-1:1
    f(j) = f(j)/L(j,j);
    for i  = max([1,j-p]):j-1
        f(i) = f(i) - L(i,j)*f(j);
    end
end
z = f;
end