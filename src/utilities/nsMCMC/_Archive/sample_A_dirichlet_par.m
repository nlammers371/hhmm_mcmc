function A_prop = sample_A_dirichlet_par(A_alpha, A_counts, n_samples)

K = size(A_alpha,1);
A_prop = NaN(K,K,n_samples);
for k = 1:K
    alpha_vec = A_alpha(:,k) + A_counts(:,k);
    % generate random column
    A_prop(:,k,:) = permute(drchrnd(alpha_vec',n_samples),[2 3 1]);
end