clear all;  close all; clc
pkg load symbolic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

limit = 0.5;           % limit drawing of arrows
BrMax = 1.0;           % remanent field [T]
DeCoeff = 0.1:0.3:1.0; % drawn multipole coefficients
Width = 1;             % Width of the cube [MagnetsNum]
ArrayRadius = 2;       % outwards radial displacement of cube
MagnetsNum = 8;        % number of cubes
Multipolarity = 1;     % multipolarity (m = 1: dipole, 2: quadrupole)
theta_grid = 0:360 / MagnetsNum:360 - 360 / MagnetsNum;
GraphsOn = 0;

if exist('tmp/') ~= 7, mkdir('tmp'); end
fn = sprintf('tmp/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

for k = 1:length(DeCoeff)

  Br = [];
  for n = 1:length(theta_grid)
    Br = [Br, antistrength(theta_grid(n), BrMax, DeCoeff(k))];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Generate an anti-array matrix with rows for cubes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate multipole coefficients of anti-array
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  zero = zeros(1, 5);

##  figure
##  hold on
##  for n = 1:5
##    quiver(0, 0, real(AntiMPoleCoeffs)(n), imag(AntiMPoleCoeffs)(n),
##    'color', MPoleColors(n),
##    'ShowArrowHead', 'off',
##    'linewidth', 3
##    )
##  end
##  hold off

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Generate a nominal array matrix with rows for cubes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate multipole coefficients of nominal array
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate multipole coefficients of degraded array
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  DegMPoleCoeffs = AntiMPoleCoeffs + NomMPoleCoeffs;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot figure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MPoleColors = ['b', 'r', 'y', 'g', 'k'];

  subplot(2, 2, k)
  hold on

  for n = 1:5
    quiver(0, 0, real(DegMPoleCoeffs)(n), imag(DegMPoleCoeffs)(n),
    'color', MPoleColors(n),
    'ShowArrowHead', 'off',
    'linewidth', 3
    )
    xlim([-0.02, 0.1])
    ylim([-0.01, 0.01])
  end

  title(sprintf('k = %d', DeCoeff(k)))
  hold off

  arg(DegMPoleCoeffs) * 360 / (2 * pi)

end

saveas(gcf, sprintf('%sspinplot', fn), 'pdf')

