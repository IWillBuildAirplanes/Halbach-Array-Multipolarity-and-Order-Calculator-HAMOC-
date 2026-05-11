% This file plots the coefficients during the degradation process

clear all;  close all; clc
pkg load symbolic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Array characteristics

BrMax = 1.0;       % remanent field [T] max
k = 0.0:0.1:1.0;  % degrading coefficient

Width = 1;         % Width of the cube
ArrayRadius = 4;   % outwards radial displacement of cube
FieldRadius = 1;   % width containing field measured (square shaped good field)
MagnetsNum = 16;    % number of cubes
Multipolarity = 1; % multipolarity (m = 1: dipole, 2: quadrupole)

%%% How many terms are calculated during multipole expansion

NomMPoleNum = 5;
AntiMPoleNum = NomMPoleNum;
DegMPoleNum = NomMPoleNum;

%%% Space grids

xGrid = -FieldRadius:0.2:FieldRadius;
yGrid = xGrid;
zGrid = 0;
ThetaGrid = 0:360 / MagnetsNum:360 - 360 / MagnetsNum;
[ReGrid, ImGrid] = meshgrid(-FieldRadius:0.2:FieldRadius);
ComplexGrid = ReGrid + ImGrid * i;
[xGridMesh, yGridMesh] = meshgrid(xGrid);
[xGridInterp, yGridInterp] = meshgrid(
  -(ArrayRadius - Width):0.05:(ArrayRadius - Width)
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
% Cycle degrade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MPoleColors = ["b", "r", "g", "m", "k"];
MPoleWidths = [3, 3, 2, 2, 1];

NomStrengths = zeros(1, MagnetsNum);

  for n = 1:MagnetsNum
    NomStrengths(n) = BrMax;
  end

NomSheets = MakeArray(
  MagnetsNum, Width, 0, NomStrengths, Multipolarity, ArrayRadius
);

NomMPoleCoeffs = MPoleElemOf(NomSheets, NomMPoleNum);

DegMPoleCoeffs = [];

for m = 1:length(k)

  AntiStrengths = AntiStrengthsOf(ThetaGrid, BrMax, k(m));

  AntiSheets = MakeArray(
    MagnetsNum, Width, 180, AntiStrengths, Multipolarity, ArrayRadius
  );

  AntiMPoleCoeffs = MPoleElemOf(AntiSheets, AntiMPoleNum);

  DegMPoleCoeffs = [DegMPoleCoeffs; AntiMPoleCoeffs + NomMPoleCoeffs];
end

DegMPoleCoeffs = DegMPoleCoeffs .* (zeros(size(DegMPoleCoeffs)) + 10) .^ 4;

DegMPoleCoeffsScaled = DegMPoleCoeffs;

for m = 1:length(k)
  DegMPoleCoeffsScaled(m, :) = DegMPoleCoeffsScaled(m, :) .* [1, 1, 10, 10, 10];
end

figure
hold on

slopes = zeros(1, 5);
for m = 1:5
  plot(k, abs(DegMPoleCoeffsScaled(:, m)), "Color", MPoleColors(m));
  slopes(m) = polyfit(k, abs(DegMPoleCoeffs(:, m)), 1)(1);
end

xlim([0, 1])
ylim([0, 400])
xlabel("k (T)")
ylabel("c_{n} G / cm^{n - 1}")

legend('Dipole', 'Quadrupole', 'Sextupole x10', 'Octopole x10', 'Decapole x10')

saveas(gcf, sprintf('%sMPoleLines', Directory), 'png')

hold off

figure
hold on

DegMPoleCoeffsArg = 360 / 2 / pi * arg(DegMPoleCoeffs);
DegMPoleCoeffsArg(1, 2) = 0;
DegMPoleCoeffsArg(1, 3) = -180;

for m = 1:5
  plot(k, DegMPoleCoeffsArg(:, m), "Color", MPoleColors(m));
end

xlim([0, 1])
ylim([-180, 180])
xlabel("k (T)")
ylabel("arg(c_{n}) (deg)")

legend('Dipole', 'Quadrupole', 'Sextupole', 'Octopole', 'Decapole')

saveas(gcf, sprintf('%sMPoleLinesArg', Directory), 'png')

hold off

