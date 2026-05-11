

function c = MPoleElemOf(sheets, num)

  indices = 0:num - 1;
  roots = exp(- 2 * pi * i * indices / num);
  rootsvec = [real(roots); imag(roots); zeros(1, num)];
  valuesvec = zeros(3, num);

  for n = indices
    valuesvec(:, n + 1) = Bsheets(sheets, rootsvec(:, n + 1));
  end

  values = valuesvec(1, :) + valuesvec(2, :) * i;

  c = ifft(- conj(i * values));
end


