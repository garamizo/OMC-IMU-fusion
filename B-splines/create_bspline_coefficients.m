function M = create_bspline_coefficients(bord)
    M = 1;
    for k = 1 : (bord-1)
        M = ([M; zeros(1,k)] * ([diag(1:(k)), zeros(k,1)] + [zeros(k,1), diag((k-1):-1:0)]) + ...
             [zeros(1,k); M] * ([-eye(k), zeros(k,1)] + [zeros(k,1), eye(k)]) ) / k;
    end
end