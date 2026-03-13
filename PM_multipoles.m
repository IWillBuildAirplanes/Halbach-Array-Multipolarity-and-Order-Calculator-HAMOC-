clear all;  close all; clc
pkg load symbolic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

limit = 0.5;       % limit drawing of arrows
BrMax = 1.0;       % remanent field [T]
k = 0.1;           % degrading coefficient
Width = 1;         % Width of the cube [MagnetsNum]
ArrayRadius = 2;   % outwards radial displacement of cube
MagnetsNum = 8;    % number of cubes
Multipolarity = 1; % multipolarity (m = 1: dipole, 2: quadrupole)
GraphsOn = 1;

if exist('tmp/') ~= 7, mkdir('tmp'); end
fn = sprintf('tmp/MagnetsNum%02d_Multipolarity%1d_w%02d_ArrayRadius%02d',
  MagnetsNum, Multipolarity, Width, ArrayRadius
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antiarray magnets' strengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_grid = 0:360 / MagnetsNum:360 - 360 / MagnetsNum;

Br = [];
for n = 1:length(theta_grid)
  Br = [Br, antistrength(theta_grid(n), BrMax, k)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate an anti-array matrix with rows for cubes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 0:MagnetsNum - 1
  cube = make_cubez(Width, 1, Br(n + 1));
  cubey = sheets_rotate_x(cube, -90);

  tmp = sheets_rotate_z(cubey, -n * Multipolarity * 360 / MagnetsNum + 180);
  tmp = sheets_translate(tmp, [ArrayRadius; 0; 0]);
  tmp = sheets_rotate_z(tmp, -n * 360 / MagnetsNum);
  if n == 0
    AntiSheets = tmp;
  else
    AntiSheets = [AntiSheets; tmp];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot anti-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GraphsOn == 1
  figure
  hold on
end

if GraphsOn == 1
  for m = 1:size(AntiSheets, 1)
    draw_sheets(AntiSheets(m, :))
  end
end

gridpoints = -(ArrayRadius - Width / 2):0.2:(ArrayRadius - Width / 2);
gx = gridpoints;  gy = gx;  gz = 0;   % in xy-plane

if GraphsOn == 1
  field_from_sheets3(gx, gy, gz, AntiSheets, limit)
end

axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate multipole coefficients of anti-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AntiMPoleCoeffs = [0, 0, 0, 0, 0];
AntiMPoleNum = length(AntiMPoleCoeffs);

BList = [];
BListVec = [0; 0; 0];
for n = 0:1:AntiMPoleNum - 1
  TestPosComp = exp(2 * pi * i * n / AntiMPoleNum);
  TestPosVec = [real(TestPosComp); imag(TestPosComp); 0];
  BListVec = Bsheets(AntiSheets, TestPosVec);
  BList = [BList, BListVec(1) + BListVec(2) * i];
end

AntiMPoleCoeffs = (fft(-conj(i * BList)) / AntiMPoleNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add field recreated from multipole coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xgrid, ygrid] = meshgrid(
  -(ArrayRadius - Width / 2):0.2:(ArrayRadius - Width / 2)
);

function B = magfield(c, z)
  B = zeros(size(z));
  for n = 1:length(c)
    B += c(n) * z .^ (n - 1);
  end
  B = -conj(B * i);
end

zgrid = xgrid + ygrid * i;
Bfield = magfield(AntiMPoleCoeffs, zgrid);
if GraphsOn == 1
  quiver(xgrid, ygrid, real(Bfield), imag(Bfield),
    'LineWidth', 2, 'Color', 'b');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a nominal array matrix with rows for cubes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 0:MagnetsNum - 1
  cube = make_cubez(Width, 1, 1);
  cubey = sheets_rotate_x(cube, -90);

  tmp = sheets_rotate_z(cubey, -n * Multipolarity * 360 / MagnetsNum);
  tmp = sheets_translate(tmp, [ArrayRadius; 0; 0]);
  tmp = sheets_rotate_z(tmp, -n * 360 / MagnetsNum);
  if n == 0
    NomSheets = tmp;
  else
    NomSheets = [NomSheets; tmp];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot nominal array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GraphsOn == 1
  figure
  hold on
end

if GraphsOn == 1
  for m = 1:size(NomSheets, 1)
    draw_sheets(NomSheets(m, :))
  end
end

gridpoints = -(ArrayRadius - Width / 2):0.2:(ArrayRadius - Width / 2);
gx = gridpoints;  gy = gx;  gz = 0;   % in xy-plane

if GraphsOn == 1
  field_from_sheets3(gx, gy, gz, NomSheets, limit)
end

axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate multipole coefficients of nominal array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NomMPoleCoeffs = [0, 0, 0, 0, 0];
NomMPoleNum = length(NomMPoleCoeffs);

BList = [];
BListVec = [0; 0; 0];
for n = 0:1:NomMPoleNum - 1
  TestPosComp = exp(2 * pi * i * n / NomMPoleNum);
  TestPosVec = [real(TestPosComp); imag(TestPosComp); 0];
  BListVec = Bsheets(NomSheets, TestPosVec);
  BList = [BList, BListVec(1) + BListVec(2) * i];
end

NomMPoleCoeffs = fft(-conj(i * BList)) / NomMPoleNum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate multipole coefficients of degraded array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DegMPoleCoeffs = AntiMPoleCoeffs + NomMPoleCoeffs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print multipole coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NomMPoleCoeffsMod = abs(NomMPoleCoeffs);
NomMPoleCoeffsArg = arg(NomMPoleCoeffs) * 360 / (2 * pi);
AntiMPoleCoeffsMod = abs(AntiMPoleCoeffs);
AntiMPoleCoeffsArg = arg(AntiMPoleCoeffs) * 360 / (2 * pi);
DegMPoleCoeffsMod = abs(DegMPoleCoeffs);
DegMPoleCoeffsArg = arg(DegMPoleCoeffs) * 360 / (2 * pi);

warning off
format short g

latex(sym(k))
latex(vpa(sym(NomMPoleCoeffsMod),3))
latex(vpa(sym(NomMPoleCoeffsArg),3))
latex(vpa(sym(AntiMPoleCoeffsMod),3))
latex(vpa(sym(AntiMPoleCoeffsArg),3))
latex(vpa(sym(DegMPoleCoeffsMod),3))
latex(vpa(sym(DegMPoleCoeffsArg),3))


