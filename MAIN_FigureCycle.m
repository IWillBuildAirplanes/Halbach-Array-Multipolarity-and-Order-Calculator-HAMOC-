% This file plots the coefficients during the degradation process

clear all;  close all; clc
pkg load symbolic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Array characteristics

BrMax = 1.0;       % remanent field [T] max
k = [0, 1];           % degrading coefficient

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

MPoleColors = ['b', 'r', 'g', 'y', 'k'];
MPoleWidths = [3, 3, 1, 1, 1];

for m = 1:2
  NomStrengths = zeros(1, MagnetsNum);

  for n = 1:MagnetsNum
    NomStrengths(n) = BrMax;
  end

  NomSheets = MakeArray(
    MagnetsNum, Width, 0, NomStrengths, Multipolarity, ArrayRadius
  );

  NomMPoleCoeffs = MPoleElemOf(NomSheets, NomMPoleNum);

  AntiStrengths = AntiStrengthsOf(ThetaGrid, BrMax, k(m));

  AntiSheets = MakeArray(
    MagnetsNum, Width, 180, AntiStrengths, Multipolarity, ArrayRadius
  );

  AntiMPoleCoeffs = MPoleElemOf(AntiSheets, AntiMPoleNum);

  DegMPoleCoeffs = AntiMPoleCoeffs + NomMPoleCoeffs;

  DegMPoleCoeffs = DegMPoleCoeffs ./ abs(DegMPoleCoeffs);

  %%% Plot

  figure
  hold on
  %title(sprintf("k = %d", k(m)))
  for n = 1:NomMPoleNum
    quiver(0, 0, real(DegMPoleCoeffs)(n), imag(DegMPoleCoeffs)(n),
      'color', MPoleColors(n),
      'ShowArrowHead', 'off',
      'LineWidth', MPoleWidths(n)
    )
  end
  xlim([-1, 1])
  ylim([-1, 1])
  axis off
  legend('Dipole', 'Quadrupole', 'Sextupole', 'Octopole', 'Decapole')
end

saveas(gcf, sprintf('%sMPoleCycle', Directory), 'png')
