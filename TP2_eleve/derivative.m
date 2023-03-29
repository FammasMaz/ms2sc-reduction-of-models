function ahat = derivative(nt, a, It)
    Ft = zeros(nt, nt);
    Fte = [-1 -1; 1 1]*0.5;
    for i = 1:nt-1
       Ft(i:i+1, i:i+1) = Ft(i:i+1, i:i+1) + Fte;
    end
    ahat = ((inv(It') * Ft')*a)';
end
