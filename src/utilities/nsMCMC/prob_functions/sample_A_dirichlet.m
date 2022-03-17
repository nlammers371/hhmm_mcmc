function A_prop = sample_A_dirichlet(A_alpha, A_counts)

K = size(A_alpha,2);
A_prop = NaN(K);
for k = 1:K
    alpha_vec = A_alpha(:,k) + A_counts(:,k);
    % generate random column
    A_prop(:,k) = drchrnd(alpha_vec',1);
end