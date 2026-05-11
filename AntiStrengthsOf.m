% antistrength.m
% Usage: gives strength of the anti (opposing) magnetic field as a function of
% magnet's angle with +x-axis which models degredation. Adding the magnet's
% nominal strength with this antistrength will be its ultimate degraded
% strength. The strenght is stronger on one side due to synchotron radiation.
% H = antistrength(theta (in degrees, BrMax (maximum strength), k (degrading
% coefficient from 0 to 1)).
function H = AntiStrengthsOf(theta, BrMax, k)
  H = zeros(1, length(theta));
  Limits = [0, 30, 90, 150, 165, 195, 210, 270, 330];
  Strengths = [0.0, 0.25, 0.5, 0.75, 1.0, 0.75, 0.5, 0.25, 0.0];
  for n = 1:length(theta)
    for m = 1:length(Limits)
      if theta(n) > Limits(m)
        theta(n);
        Limits(m);
        H(n) = Strengths(m) * k * BrMax;
      end
  end
end
