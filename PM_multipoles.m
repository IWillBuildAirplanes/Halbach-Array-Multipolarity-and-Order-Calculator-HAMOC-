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
xGrid = -(ArrayRadius - Width / 2):0.2:(ArrayRadius - Width / 2);
yGrid = xGrid;
zGrid = 0;

if exist('tmp/') ~= 7, mkdir('tmp'); end
fn = sprintf('tmp/MagnetsNum%02d_Multipolarity%1d_w%02d_ArrayRadius%02d',
  MagnetsNum, Multipolarity, Width, ArrayRadius
);


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
  title("Nominal Array Field")

  hold on

  for m = 1:size(NomSheets, 1)
    draw_sheets(NomSheets(m, :))
  end

  [NomFieldX, NomFieldY, NomFieldZ] = field_from_sheets3(xGrid, yGrid, zGrid, NomSheets, limit, 'r');

  axis equal
  view(0, 90)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate multipole coefficients of nominal array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NomMPoleNum = 5;
NomMPoleCoeffs = zeros(1, NomMPoleNum);

NomFieldTestValues = [];
NomFieldTestEntry = zeros(3, 1);

for n = 0:1:NomMPoleNum - 1
  TestPosComplex = exp(2 * pi * i * n / NomMPoleNum);
  TestPosVector = [real(TestPosComplex); imag(TestPosComplex); 0];
  NomFieldTestEntry = Bsheets(NomSheets, TestPosVector);
  NomFieldTestValues = [
    NomFieldTestValues,
    NomFieldTestEntry(1) + NomFieldTestEntry(2) * i
  ];
end

NomMPoleCoeffs = fft(-conj(i * NomFieldTestValues)) / NomMPoleNum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antiarray magnets' strengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ThetaGrid = 0:360 / MagnetsNum:360 - 360 / MagnetsNum;

Br = [];
for n = 1:length(ThetaGrid)
  Br = [Br, antistrength(ThetaGrid(n), BrMax, k)];
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
  title("Antiarray Field (Red) and Recreated Field from MPole Coefficients (Blue)")

  hold on

  for m = 1:size(AntiSheets, 1)
    draw_sheets(AntiSheets(m, :))
  end

  AntiField = field_from_sheets3(xGrid, yGrid, zGrid, AntiSheets, limit, 'r');
  axis equal
  view(0, 90)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate multipole coefficients of anti-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AntiMPoleNum = 5;
AntiMPoleCoeffs = zeros(1, AntiMPoleNum);

AntiFieldTestValues = [];
AntiFieldTestEntry = zeros(3, 1);

for n = 0:1:AntiMPoleNum - 1
  TestPosComplex = exp(2 * pi * i * n / AntiMPoleNum);
  TestPosVector = [real(TestPosComplex); imag(TestPosComplex); 0];
  AntiFieldTestEntry = Bsheets(AntiSheets, TestPosVector);
  AntiFieldTestValues = [
    AntiFieldTestValues,
    AntiFieldTestEntry(1) + AntiFieldTestEntry(2) * i
  ];
end

AntiMPoleCoeffs = (fft(-conj(i * AntiFieldTestValues)) / AntiMPoleNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add antiarray field recreated from multipole coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = FieldFromCoeffs(c, z)
  B = zeros(size(z));
  for n = 1:length(c)
    B += c(n) * z .^ (n - 1);
  end
  B = -conj(B * i);
end

[ReGrid, ImGrid] = meshgrid(
  -(ArrayRadius - Width / 2):0.2:(ArrayRadius - Width / 2)
);

ComplexGrid = ReGrid + ImGrid * i;
RecreatedField = FieldFromCoeffs(AntiMPoleCoeffs, ComplexGrid);

if GraphsOn == 1
  quiver(xGrid, yGrid, real(RecreatedField), imag(RecreatedField),
    'LineWidth', 2, 'color', 'b');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate multipole coefficients of degraded array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DegMPoleCoeffs = AntiMPoleCoeffs + NomMPoleCoeffs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot degraded array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GraphsOn == 1
  figure
  title("Nominal Array (Red) versus Degraded Array (Black)")

  hold on

  for m = 1:size(NomSheets, 1)
    draw_sheets(NomSheets(m, :))
  end

  field_from_sheets3(xGrid, yGrid, zGrid, NomSheets, limit, 'r');
  axis equal
  view(0, 90)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add degraded field recreated from multipole coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DegFieldRecreated = FieldFromCoeffs(DegMPoleCoeffs, ComplexGrid);

if GraphsOn == 1
  quiver(xGrid, yGrid, real(DegFieldRecreated), imag(DegFieldRecreated),
    'LineWidth', 2, 'color', 'k');
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

GFRIndices = [];
for i = 1:length(xGrid)
  if abs(xGrid(i)) <= 1
    GFRIndices = [GFRIndices, i];
  end
end

NomField = NomFieldX + NomFieldY * i;
Error = abs(NomField - DegFieldRecreated);
xGridInterp = -(ArrayRadius - Width / 2):0.05:(ArrayRadius - Width / 2);
yGridInterp = xGridInterp;
% ErrorInterp = interp2(xGrid, yGrid, Error);
surf(xGrid, yGrid, Error)

