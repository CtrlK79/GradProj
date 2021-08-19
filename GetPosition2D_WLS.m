function p = GetPosition2D_WLS(Planets, EMatrix, p_prev)
    n = length(Planets(1, :));
    E = zeros(2, 2, n);
    Er = zeros(2, n);
    
    for i = 1:n
        E(:, :, i) = (eye(2) - EMatrix(:, i)*EMatrix(:, i)')/(norm(Planets(:, i)-p_prev))^2;
        Er(:, i) = E(:, :, i) * Planets(:, i);
    end
    
    p = sum(E, 3)\sum(Er, 2);
end