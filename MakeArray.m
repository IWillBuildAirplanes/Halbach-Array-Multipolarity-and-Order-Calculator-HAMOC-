

function array = MakeArray(Num, Width, Turn, Strengths, Multipolarity, Radius)
  array = [];
  for n = 0:Num - 1
    cube = make_cubez(Width, 1, Strengths(n + 1));
    cubey = sheets_rotate_x(cube, -90);
    tmp = sheets_rotate_z(cubey, -n * Multipolarity * 360 / Num + Turn);
    tmp = sheets_translate(tmp, [Radius; 0; 0]);
    tmp = sheets_rotate_z(tmp, -n * 360 / Num);
    array = [array; tmp];
  end
end

