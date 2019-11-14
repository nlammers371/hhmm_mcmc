function A_prop = sample_A(A_curr,n)
    A_prop = zeros(size(A_curr));
    for k = 1:size(A_curr,1)
        A_prop(:,k) = mnrnd(n,A_curr(:,k));
    end
    A_prop = A_prop ./ sum(A_prop);
        