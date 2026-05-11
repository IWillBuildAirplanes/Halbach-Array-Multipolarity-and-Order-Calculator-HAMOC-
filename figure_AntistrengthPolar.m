thetagrid = 0:360;
antistrengthgrid = [];
if exist('figures/') ~= 7
  mkdir('figures');
end

Directory = sprintf('figures/');
for n = thetagrid
  antistrengthgrid = [antistrengthgrid, AntiStrengthsOf(n, 1, 1)];
end
figure
polar(thetagrid * pi / 180, antistrengthgrid)
saveas(gcf, sprintf('%sAntiStrength', Directory), 'png')
