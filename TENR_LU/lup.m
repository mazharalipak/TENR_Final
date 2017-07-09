function [L,U,P] = lup(A)
% lup factorization with partial pivoting
% [L,U,P] = lup(A) returns unit lower triangular matrix L, upper
% triangular matrix U, and permutation matrix P so that P*A = L*U.
    n = length(A);
    L = zeros(n);
    U = zeros(n);
    P = eye(n);


    for k=1:n
        % find the entry in the left column with the largest abs value (pivot)
        [~,r] = max(abs(A(k:end,k)));
        r = n-(n-k+1)+r;    

        A([k r],:) = A([r k],:);
        P([k r],:) = P([r k],:);
        L([k r],:) = L([r k],:);

        % from the pivot down divide by the pivot
        L(k:n,k) = A(k:n,k) / A(k,k);

        U(k,1:n) = A(k,1:n);
        A(k+1:n,1:n) = A(k+1:n,1:n) - L(k+1:n,k)*A(k,1:n);

    end
    U(:,end) = A(:,end);

end