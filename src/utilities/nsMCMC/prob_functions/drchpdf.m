function prob = drchpdf(x,a)


% take a sample from a dirichlet distribution
numerator = prod(x.^(a-1));
prefactor = gamma(sum(a))/prod(gamma(a));
prob = prefactor*numerator;