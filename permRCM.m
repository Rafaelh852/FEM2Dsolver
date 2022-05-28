function [AA,Q] = permRCM(A)
%A is a square
% AA is permuted matrix QAQ' 

p = symrcm(A); % permutation vector p

AA= A(p,p); %reorders A to have diagonal elements closer to diagonal than A. 

% need to make p a matrix. ( we call it Q)
Q = zeros(length(A(:,1)), length(A(:,1)));

% fill matrix Q with ones in the p index. 
for row = 1:length(p)
    Q(row,p(row)) = 1;
end


Q = sparse(Q);
end