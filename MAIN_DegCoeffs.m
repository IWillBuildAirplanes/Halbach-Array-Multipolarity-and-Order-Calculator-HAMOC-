% This file gets multipole coefficients of degraded arrays

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Degraded array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate multipole coefficients of degraded array

% Multipole expansion of magnetic field is linear in multipole coefficients

DegMPoleCoeffs = AntiMPoleCoeffs + NomMPoleCoeffs;

%%% Plot degraded field recreated from multipole coefficients

DegFieldRecreated = FieldFromCoeffs(DegMPoleCoeffs, ComplexGrid);

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
format shortG

latex(sym(k))
display("Nominal Coefficients")
latex(vpa(sym(NomMPoleCoeffsMod * 10000), 3))
latex(vpa(sym(NomMPoleCoeffsArg), 3))
display("Degraded Coefficients")
latex(vpa(sym(DegMPoleCoeffsMod * 10000) ,3))
latex(vpa(sym(DegMPoleCoeffsArg), 3))
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##% Plot error
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
##if GraphsOn == 1
##  figure
##  title("Error between Nominal and Degraded Field (T)")
##  xlabel("x-axis (cm)")
##  ylabel("y-axis (cm)")
##  hold on
##
##  view(0, 90)
##  xlim([-1, 1])
##  ylim([-1, 1])
##end
##
##NomFieldMagInterp = interp2(
##  xGridMesh, yGridMesh, abs(NomField), xGridInterp, yGridInterp
##);
##
##DegFieldMagInterp = interp2(
##  xGridMesh, yGridMesh, abs(DegField), xGridInterp, yGridInterp
##);
##
##FieldError(ComplexGridInterp, NomFieldMagInterp, DegFieldMagInterp)
##
##saveas(gcf, sprintf('%sDegradedArrayFieldError', Directory), 'png')
