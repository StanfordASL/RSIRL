function w_new = tighten(w, A, b, L)

P = Polyhedron('V', eye(L)) & Polyhedron('A', A, 'b', b-w);
P = P.minVRep();
V = P.V; 
% n_V = size(V,1);
M = size(A,1);
b_ = zeros(M,1);

for i=1:M
%     sat_vertices = 0;
    b_(i,1) = max(V*A(i,:)');
%     for j=1:n_V
%         if abs(V(j,:)*A(i,:)' - b_(i,1)) <=0.0001
%             sat_vertices = sat_vertices + 1;
%         end
%     end
%     if sat_vertices == 1
%         b_(i,1) = b_(i,1) - 0.01;
%     end
    if (b(i)-w(i) > b_(i) + 0.001)
        b_(i) = b_(i) - 0.001;
    end
end
w_new = b - b_;
    
end


