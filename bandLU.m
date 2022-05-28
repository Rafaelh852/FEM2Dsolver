function [K] = bandLU(K,p,q)

% test K
%K= [1 2 0 0;3 1 2 2;0 0 1 3;0 0 1 1];
%[p,q] = bandwidth(K);

n = length(K);

for m = 1:n-1
    for i = m+1:min([m+p,n])
        K(i,m) = K(i,m)/K(m,m);
    end
    for j = m+1:min([m+q,n])
        for i = m+1:min([m+p,n])
            K(i,j)= K(i,j) - K(i,m)*K(m,j);
        end
    end
end

K = K;
% %original K == (tril(K,-1) + eye(n))*triu(K)
% K
% L = tril(K,-1) + eye(n)
% U = triu(K)

% %test vector
% f = [2; 2; 3; 4];
% 
% %--- Bandforward ---
% for j = 1:n
%     for i= j+1:min([j+p,n])
%         f(i) = f(i) - L(i,j)*f(j);
%     end 
% end
% 
% %checking solution
% f
% L*f
% 
% %--- Bandbackward ---    
% for j = n:-1:1
%     f(j) = f(j)/U(j,j);
%     for i  = max([1,j-q]):j-1
%         f(i) = f(i) - U(i,j)*f(j);
%     end
% end
% 
% %checking solution
% f
% U*f

end