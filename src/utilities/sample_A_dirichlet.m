function A_prop = sample_A_dirichlet(A_alpha, A_counts,n_reps)

    K = size(A_alpha,1);
    A_prop = NaN(K,K,n_reps);
    for k = 1:K
        alpha_vec = A_alpha(:,k) + A_counts(:,k);
        % generate random column
        A_prop(:,k,:) = reshape(drchrnd(alpha_vec',n_reps)',K,1,[]);
    end