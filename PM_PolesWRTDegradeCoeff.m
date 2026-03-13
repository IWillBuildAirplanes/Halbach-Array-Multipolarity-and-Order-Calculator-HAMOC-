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
Directory = sprintf('tmp/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate multipole values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize k (Degerading Coefficient) and Multipole Element Lists

DeCoeff = 0:0.05:1.0;

Dipole = [];
Quadrupole = [];
Sextupole = [];
Octopole = [];
Decapole = [];

DipoleAngle = [];
QuadrupoleAngle = [];
SextupoleAngle = [];
OctopoleAngle = [];
DecapoleAngle = [];

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

  Dipole = [Dipole, DegMPoleCoeffs(1)];
  Quadrupole = [Quadrupole, DegMPoleCoeffs(2)];
  Sextupole = [Sextupole, DegMPoleCoeffs(3)];
  Octopole = [Octopole, DegMPoleCoeffs(4)];
  Decapole = [Decapole, DegMPoleCoeffs(5)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All element plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

hold on

plot(DeCoeff, abs(Dipole),
  'color', 'b',
  'linewidth', 3,
  'DisplayName', 'Dipole'
);
plot(DeCoeff, abs(Quadrupole),
  'color', 'r',
  'linewidth', 1.5,
  'DisplayName', 'Quadrupole'
);
plot(DeCoeff, abs(Sextupole),
  'color', 'y',
  'linewidth', 1.5,
  'DisplayName', 'Sextupole'
);
plot(DeCoeff, abs(Octopole),
  'color', 'g',
  'linewidth', 1.5,
  'DisplayName', 'Octopole'
);
plot(DeCoeff, abs(Decapole),
  'color', 'k',
  'linewidth', 1.5,
  'DisplayName', 'Decapole'
);

title("Multiple Elements' Magnitude wrt Degrading Coefficient")

legend_names = {'Dipole', 'Quadrupole', 'Sextupole', 'Octopole', 'Decapole'};

legend(legend_names, 'location', 'northwest')

saveas(gcf, sprintf('%scoolplot', Directory), 'pdf')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoomed in plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

hold on

plot(DeCoeff, Quadrupole,
  'color', 'r',
  'linewidth', 2
);
plot(DeCoeff, Sextupole,
  'color', 'y',
  'linewidth', 2
);
plot(DeCoeff, Octopole,
  'color', 'g',
  'linewidth', 2
);
plot(DeCoeff, Decapole,
  'color', 'k',
  'linewidth', 2
);

title("Multiple Elements' Magnitude (no Dipole) wrt Degrading Coefficient")

legend_names = {'Quadrupole ', 'Sextupole', 'Octopole', 'Decapole'};

legend(legend_names, 'location', 'northwest')

saveas(gcf, sprintf('%scoolzoomplot', Directory), 'pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Angles plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

hold on

plot(DeCoeff, arg(Dipole),
  'color', 'b',
  'linewidth', 2
);
plot(DeCoeff, arg(Quadrupole),
  'color', 'r',
  'linewidth', 2
);
plot(DeCoeff, arg(Sextupole),
  'color', 'y',
  'linewidth', 2
);
plot(DeCoeff, arg(Octopole),
  'color', 'g',
  'linewidth', 2
);
plot(DeCoeff, arg(Decapole),
  'color', 'k',
  'linewidth', 2
);

% Linear regression slopes

m1 = polyfit(DeCoeff(5:21), Dipole(5:21), 1)
m2 = polyfit(DeCoeff(5:21), Quadrupole(5:21), 1)
m3 = polyfit(DeCoeff(5:21), Sextupole(5:21), 1)
m4 = polyfit(DeCoeff(5:21), Octopole(5:21), 1)
m5 = polyfit(DeCoeff(5:21), Decapole(5:21), 1)

title("Multiple Elements' Angles wrt Degrading Coefficient")

legend_names = {
  sprintf('Dipole, m = %.3f', m1(2)),
  sprintf('Quadrupole, m = %.3f', m2(2)),
  sprintf('Sextupole, m = %.3f', m3(2)),
  sprintf('Octopole, m = %.3f', m4(2)),
  sprintf('Decapole, m = %.3f', m5(2))
};

legend(legend_names, 'location', 'northwest')

saveas(gcf, sprintf('%sangleplot', Directory), 'pdf')

