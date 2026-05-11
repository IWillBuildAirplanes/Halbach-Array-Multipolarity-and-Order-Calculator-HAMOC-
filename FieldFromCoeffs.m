% fFieldFromCoeffs.m
% Usage: Calculates field at complex coordinate (z) according to multipole
% coefficients c
function B = FieldFromCoeffs(c, z)
  C = zeros(size(z));
  for n = 1:length(c)
    C += c(n) * z .^ (n - 1);
  end
  B = -conj(C * i);
end
