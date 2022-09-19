function logprob = logdrchpdf(x,a)


% take a sample from a dirichlet distribution
numerator = sum((a-1).*log(x),2);
prefactor = gammaln(sum(a,2))-sum(gammaln(a),2);
logprob = prefactor+numerator;