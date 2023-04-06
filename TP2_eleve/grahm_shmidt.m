function Q = grahm_shmidt(A, K)
% A: an mxn matrix
% Q: an mxn matrix with orthonormal columns
[m, n] = size(A);
Q = zeros(m, n);

for j = 1:n
    v = A(:, j);
    for i = 1:j-1
        q = Q(:, i);
        v = v - q*(q'*v);
    end
    Q(:, j) = v/norm(v);
end
Q = Q/sqrtm(Q'*K*Q);
end