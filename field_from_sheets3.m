% field_from_sheets3.m
% Usage: Plots the field from a set of magnetic sheets across the grid
% 'gx', 'gy', 'gz' with vectors of length scaled by 'length' and color 'color'.
function [BX, BY, BZ] = field_from_sheets3(gx, gy, gz, sheets, thick, len, col)

[XX, YY, ZZ] = meshgrid(gx, gy, gz);

BX = 0 * XX;
BY = 0 * YY;
BZ = 0 * ZZ;

%n=length(x);
for kx = 1:length(gx)
  for ky = 1:length(gy)
    for kz = 1:length(gz)
      xx = XX(kx, ky, kz);
      yy = YY(kx, ky, kz);
      zz = ZZ(kx, ky, kz);
      %if (abs(xx)>a2) | (abs(yy)>a2) | (abs(zz)>a2)
        r2 = [xx; yy; zz];
        B = Bsheets(sheets, r2);
        BX(kx, ky, kz) = B(1);
        BY(kx, ky, kz) = B(2);
        BZ(kx, ky, kz) = B(3);
      %end
    end
  end
end

quiver3(XX, YY, ZZ, BX, BY, BZ, len, 'LineWidth', thick, 'Color', col)

xlabel('x');
ylabel('y');
zlabel('z');

set(gca,'FontSize',16)
