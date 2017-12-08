%halfplanes
%L = L+1; %virtual continuation mode
% A = zeros(L);
I_L = eye(L);
% for j = 1:L
%     n = (I_L(:,j) - (1/L)*ones(L,1));
%     A(j,:) = (n/norm(n))';
% end
A = I_L;
A = [A; -A];

M = size(A,1);

% linear parameterization: A*v <= b - w
% compute offsets when w = 0 (minimal offsets such that proba simplex included in parameterized polytope)
I_M = eye(M);
for i=1:M
    a = A(i,:)';
    V = eye(L);
    res = V*a;
    b(i,1) = max(res);
end

%% Define polytope of feasible parameters
%A*q + w <= b
% -q <= 0
%- w <= -0.02
%sum(q) == 1

A_ = [A, I_M; 
      - I_L, zeros(L,M); 
      -zeros(M,L), -I_M];
  
b_ = [b; zeros(L,1); - 0.02*ones(M,1)]; A_e = [ones(1,L), zeros(1,M)]; b_e = [1.];
P_proj = Polyhedron('A', A_, 'b', b_, 'Ae', A_e, 'be', b_e);
P_proj = P_proj.projection(L+1:L+M);

%setup LP solver
%LP_prob = setup_LP_tomlab(L,A,b);
LP_prob = setup_LP_yalmip(L,A,b);

%call once to setup warm start
LP_sol = struct('sol',0);
%[~,~,LP_sol] = compute_LP_tomlab(LP_prob,randn(L,1),kron(ones(1,M),b),M,L,LP_sol);
[~,~,LP_sol] = compute_LP_yalmip(LP_prob,randn(L,1),kron(ones(1,M),b),M,L);


