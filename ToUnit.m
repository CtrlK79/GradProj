function U = ToUnit(A)
    n = length(A(1, :));
    for i = 1:n
        U(:, i) = A(:, i)/norm(A(:, i));
    end
end