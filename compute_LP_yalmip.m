function [rho, q_opt_mat, sol] = compute_LP_yalmip(Prob,costs,b_U,M,L)
%Inputs:
%costs: LP cost vector
%b_U: parameterized halfplane offsets for each perturbation
%warm: previous tomlab solution struct

b_global = reshape(b_U,M*M,1);
costs_global = kron(ones(M,1),costs);
% try
[q_opt_vec,~] = Prob(b_global,-costs_global);
% catch
%     keyboard;
% end
q_opt_mat = reshape(q_opt_vec,L,M);

% if (~infeas)
sol.sol = 1;
rho = q_opt_mat'*costs;

end
