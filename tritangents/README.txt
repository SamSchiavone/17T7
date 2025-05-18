Step 1: in Julia, run `myhomotopy.jl` to find 720 smooth solutions using 
numerical homotopy, refine solutions using `refine.jl` to desired precision
(30000 is overkill), save output to `sols_refined.out`.  Complete run logged in `homotopy.log`.

Step 2: run `regex-m-jl.py` in python or sage to convert to magma, output `sols_refined.m`.

Step 3: run `exact-pols.m` to produce exact candidate equations for the
coefficients of the plane (degree 120 polynomials) and verify that 
the plane is indeed tritangent by an exact calculation with resultants.