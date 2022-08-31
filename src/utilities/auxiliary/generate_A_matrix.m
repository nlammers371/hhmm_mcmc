function A = generate_A_matrix(Q)

  diag_terms = diag(-Q);
  p_jump = 1-exp(-diag_terms);
  A = zeros(3);
  A(eye(3)==1) = 1-p_jump;
  A(2,1) = 1-A(1,1);
  A(2,3) = 1-A(3,3);
  pj2 = 1-A(2,2);
  k_rat = Q(1,2)/(Q(3,2)+Q(1,2));
  A(1,2) = pj2*k_rat;
  A(3,2) = pj2*(1-k_rat);