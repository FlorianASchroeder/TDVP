function [U, S, V] = svd2(T)
% SVD and returns U, S, V' for any matrix.
% Transposes T in case m < n:  to make code shorter
%

[m, n] = size(T);
if m >= n
    [U, S, V] = svd(T, 0);
else
    [V, S, U] = svd(T', 0);
end

V = V';