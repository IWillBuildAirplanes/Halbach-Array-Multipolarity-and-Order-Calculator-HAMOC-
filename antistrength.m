% antistrength.m
% Usage: gives strength of the anti (opposing) magnetic field as a function of
% magnet's angle with +x-axis which models degredation. Adding the magnet's
% nominal strength with this antistrength will be its ultimate degraded
% strength. The strenght is stronger on one side due to synchotron radiation.
% H = antistrength(theta (in degrees, BrMax (maximum strength), k (degrading
% coefficient from 0 to 1)).
function H = antistrength(theta, BrMax, k)
  if theta < 30
  H = 0.0;
  elseif theta < 90
  H = 0.25 * k * BrMax;
  elseif theta < 150
  H = 0.5 * k * BrMax;
  elseif theta < 165
  H = 0.75 * k * BrMax;
  elseif theta < 180
  H = 1.0 * k * BrMax;
  elseif theta < 195
  H = 1.0 * k * BrMax;
  elseif theta < 210
  H = 0.75 * k * BrMax;
  elseif theta < 270
  H = 0.5 * k * BrMax;
  elseif theta < 330
  H = 0.25 * k * BrMax;
  else
  H = 0;
  end
end
