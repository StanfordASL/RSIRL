function w_tight = w_tightest(A,b,P)

V = P.V;
M = size(A,1);
for i=1:M
    w_tight(i,1) = b(i,1) - max(V * A(i,:)');
end

end