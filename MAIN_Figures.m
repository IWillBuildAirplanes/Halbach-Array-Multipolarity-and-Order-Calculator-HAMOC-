% This file creates a whole variety of figures related to the Halbach array
% degradation. Basically, it creates a nominal array and an anti-array, the
% superposition of which creates the degraded array.

clear all;  close all; clc
pkg load symbolic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Array characteristics

BrMax = 1.0;       % remanent field [T] max
k = 1.0;           % degrading coefficient

Width = 1;         % Width of the cube
ArrayRadius = 4;   % outwards radial displacement of cube
FieldRadius = 1;   % width containing field measured (square shaped good field)
MagnetsNum = 16;    % number of cubes
Multipolarity = 1; % multipolarity (m = 1: dipole, 2: quadrupole)

%%% How many terms are calculated during multipole expansion

NomMPoleNum = 5;
AntiMPoleNum = NomMPoleNum;
DegMPoleNum = NomMPoleNum;

%%% Display stuff

VecMagThick = 2.0; % line thickness of vectors
VecMagScale = 0.8; % visual scaling of vectors in plots

%%% Space grids

xGrid = -FieldRadius:0.2:FieldRadius;
yGrid = xGrid;
zGrid = 0;
ThetaGrid = 0:360 / MagnetsNum:360 - 360 / MagnetsNum;
[ReGrid, ImGrid] = meshgrid(-FieldRadius:0.2:FieldRadius);
ComplexGrid = ReGrid + ImGrid * i;
[xGridMesh, yGridMesh] = meshgrid(xGrid);
[xGridInterp, yGridInterp] = meshgrid(
  -FieldRadius:0.05:FieldRadius
);
ComplexGridInterp = interp2(
  xGrid, yGrid, ComplexGrid, xGridInterp, yGridInterp
);

%%% Directory

if exist('figures/') ~= 7
  mkdir('figures');
end

Directory = sprintf('figures/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nominal Array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Generate array

NomStrengths = zeros(1, MagnetsNum);

for n = 1:MagnetsNum
  NomStrengths(n) = BrMax;
end

NomSheets = MakeArray(
  MagnetsNum, Width, 0, NomStrengths, Multipolarity, ArrayRadius
);

%%% Plot array

figure
hold on

% Plot each sheet of each cube of array

for m = 1:size(NomSheets, 1)
  draw_sheets(NomSheets(m, :))
end

% Plot field of nominal array

[NomFieldX, NomFieldY, NomFieldZ] = field_from_sheets3(
  xGrid, yGrid, zGrid,
  NomSheets,
  VecMagThick, VecMagScale,
  'b'
);

NomField = NomFieldX + NomFieldY * i;

axis off
view(20, 75)
xlim([-5, 5])
ylim([-5, 5])
zlim([-1, 1])
%saveas(gcf, sprintf('%sNominalArrayFieldWhole', Directory), 'png')

axis on
view(0, 90)
xlim([-1.2, 1.2])
ylim([-1.2, 1.6])
xlabel("x-axis (cm)")
ylabel("y-axis (cm)")
saveas(gcf, sprintf('%sNominalArrayField', Directory), 'png')

%%% Calculate multipole coefficients

NomMPoleCoeffs = MPoleElemOf(NomSheets, NomMPoleNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anti-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Strength of each magnet of antiarray as function of azimuthal angle

AntiStrengths = AntiStrengthsOf(ThetaGrid, BrMax, k);

%%% Generate anti-array

AntiSheets = MakeArray(
  MagnetsNum, Width, 180, AntiStrengths, Multipolarity, ArrayRadius
);

%%% Calculate multipole coefficients of anti-array

AntiMPoleCoeffs = MPoleElemOf(AntiSheets, AntiMPoleNum);

%%% Plot anti-array

figure
hold on

% Plot field of anti-array

[AntiFieldX, AntiFieldY, AntiFieldZ] = field_from_sheets3(
  xGrid, yGrid, zGrid,
  AntiSheets,
  3, VecMagScale,
  'r'
);

AntiField = AntiFieldX + AntiFieldY * i;

% Recreate field from coefficients

AntiFieldRecreated = FieldFromCoeffs(AntiMPoleCoeffs, ComplexGrid);

quiver(xGrid, yGrid, real(AntiFieldRecreated), imag(AntiFieldRecreated),
  VecMagScale, 'LineWidth', 1, 'color', '#00E020');

axis equal
view(0, 90)
xlim([-1.2, 1.2])
ylim([-1.2, 1.2])
xlabel("x-axis (cm)")
ylabel("y-axis (cm)")
saveas(gcf, sprintf('%sAntiArrayFieldRecreated', Directory), 'png')

%%% Recreation error

AntiFieldInterp = interp2(
  xGridMesh, yGridMesh, AntiField, xGridInterp, yGridInterp
);

AntiFieldRecreatedInterp = interp2(
  xGridMesh, yGridMesh, AntiFieldRecreated, xGridInterp, yGridInterp
);

FieldErrorMag(ComplexGridInterp, AntiFieldInterp, AntiFieldRecreatedInterp)

view(0, 90)
xlim([-1, 1])
ylim([-1. 1])
xlabel("x-axis (cm)")
ylabel("y-axis (cm)")
saveas(gcf, sprintf('%sAntiArrayFieldRecreatedErrorMag', Directory), 'png')

FieldErrorArg(ComplexGridInterp, AntiFieldInterp, AntiFieldRecreatedInterp)

view(0, 90)
xlim([-1, 1])
ylim([-1. 1])
xlabel("x-axis (cm)")
ylabel("y-axis (cm)")
saveas(gcf, sprintf('%sAntiArrayFieldRecreatedErrorArg', Directory), 'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Degraded array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate multipole coefficients of degraded array

% Multipole expansion of magnetic field is linear in multipole coefficients

DegMPoleCoeffs = AntiMPoleCoeffs + NomMPoleCoeffs;

%%% Plot nominal versus degraded field

figure
hold on

% Nominal field

field_from_sheets3(
  xGrid, yGrid, zGrid,
  NomSheets,
  VecMagThick, VecMagScale,
  'r'
);

% Plot degraded field recreated from multipole coefficients

DegFieldRecreated = FieldFromCoeffs(DegMPoleCoeffs, ComplexGrid);

##DegMPoleCoeffs = MPoleElemOf([NomSheets; AntiSheets], DegMPoleNum)

quiver(xGrid, yGrid,
  real(DegFieldRecreated), imag(DegFieldRecreated),
  VecMagScale, 'LineWidth', 1.5, 'color', 'k'
);

axis equal
view(0, 90)
xlim([-1.2, 1.2])
ylim([-1.2, 1.6])
xlabel("x-axis (cm)")
ylabel("y-axis (cm)")
saveas(gcf, sprintf('%sDegradedArrayFieldMax', Directory), 'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NomFieldInterp = interp2(
  xGridMesh, yGridMesh, NomField, xGridInterp, yGridInterp
);

DegFieldInterp = interp2(
  xGridMesh, yGridMesh, DegFieldRecreated, xGridInterp, yGridInterp
);

FieldErrorMag(ComplexGridInterp, NomFieldInterp, DegFieldInterp)

view(0, 90)
xlim([-1, 1])
ylim([-1, 1])
xlabel("x-axis (cm)")
ylabel("y-axis (cm)")
saveas(gcf, sprintf('%sDegradedArrayFieldErrorMag', Directory), 'png')

FieldErrorArg(ComplexGridInterp, NomFieldInterp, DegFieldInterp)

view(0, 90)
xlim([-1, 1])
ylim([-1, 1])
xlabel("x-axis (cm)")
ylabel("y-axis (cm)")
saveas(gcf, sprintf('%sDegradedArrayFieldErrorArg', Directory), 'png')
