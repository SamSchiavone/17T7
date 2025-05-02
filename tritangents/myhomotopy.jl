# from https://www.juliahomotopycontinuation.org/examples/tritangents/ #

using HomotopyContinuation, LinearAlgebra
@var h[1:3]
@var x[1:3] y[1:3] z[1:3]
setprecision(1024)
Q = -8*x[1]^2 + 8*x[1]*x[2] - 34*x[1]*x[3] - 10*x[1] + 17*x[2]^2 - 2*x[2]*x[3] -
          9*x[2] - 28*x[3]^2 - 18*x[3] + 2
C = 4*x[1]^3 - 6*x[1]^2*x[2] + 12*x[1]^2*x[3] + 2*x[1]^2 - 6*x[1]*x[2]^2 +
          6*x[1]*x[2]*x[3] + 7*x[1]*x[2] - 12*x[1]*x[3]^2 + 4*x[1]*x[3] - 20*x[1] +
          24*x[2]^2*x[3] + 4*x[2]^2 - 13*x[2]*x[3] - 24*x[3]^3 - 8*x[3]^2 - 3*x[3] -
          12
P_x = [
        h ⋅ x - 1;
        Q;
        C;
        det([h differentiate(Q, x) differentiate(C, x)])
      ];
P_y = [p([h; x] => [h; y]) for p in P_x];
P_z = [p([h; x] => [h; z]) for p in P_x];
F = System([P_x; P_y; P_z]; variables = [h;x;y;z])
S = solve(F)
