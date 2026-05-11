

function FieldErrorMag(grid, field1, field2)
  dx = grid(1, 2) - grid(1, 1);
  dy = abs(grid(2, 1) - grid(1, 1));
  field1re = real(field1);
  field1im = imag(field1);
  field1mag = abs(field1);
  field1arg = arg(field1);

  field2re = real(field2);
  field2im = imag(field2);
  field2mag = abs(field2);
  field2arg = arg(field2);

  figure
  errmag = 100 * (field1mag - field2mag) ./ field1mag;
  surf(real(grid), imag(grid), errmag)
  cb = colorbar;
  title(cb, '## %')
  errortot = 0.0;
  for i = 1:length(errmag(1, :))
    for j = 1:length(errmag(:, 1))
      errortot += errmag(j, i) * dx * dy;
    end
  end
  errmagavg = errortot / 4
end
