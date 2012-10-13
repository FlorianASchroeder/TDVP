function [U, S, V] = svd2(T)

[m, n] = size(T);
if m >= n
    [U, S, V] = svd(T, 0);
else
    [V, S, U] = svd(T', 0);
end

V = V';