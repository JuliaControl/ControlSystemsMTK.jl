using ControlSystems, ControlSystemsMTK, ModelingToolkit
function tsmmv(J0, J1, k0, k1, c0, c1, c2)
    A = [
        0 1 0 0 0 0
        0 -c0 0 0 0 0 
        0 0 0 1 0 0
        k0 0 -(k1+k0) -c1 k1 c2
        0 0 0 0 0 1
        0 0 k1 0 -k1 -c2
    ] ./ [1, 1, 1, J0, 1, J1]
    B = [0, 1, 0, 0, 0, 0]
    C = [
        k0 0 -(k1+k0) -c1 k1 c2
    ]
    ss(A,B,C,0)
end

J0,   J1,      k0,  k1,     c0,   c1,    c2,  g =
14.5, 0.25, 28012, 4556,  2307,  155, 0.005, 103
syst1 = g*tsmmv(J0, J1, k0, k1, c0, c1, c2)
syst2 = g*tsmmv(J0, J1+2.3, k0, k1, c0, c1, c2)
syst3 = g*tsmmv(J0, J1+7.5, k0, k1, c0, c1, c2)

@variables Jm J0 J1 k0 k1 c0 c1 c2 s g bf
vars = [J0, J1,  k0,  k1, c0, c1,  c2, g]
x0   = [14.5,  0.25, 28012, 4556,  2307,  155, 0.005, 103] 

sys1 = g*tsmmv(J0, J1, k0, k1, c0, c1, c2)
sys2 = g*tsmmv(J0, J1+2.3, k0, k1, c0, c1, c2)
sys3 = g*tsmmv(J0, J1+7.5, k0, k1, c0, c1, c2)

freqs = 2pi .* exp10.(LinRange(-2, log10(100), 15))
cost, funs = ControlSystemsMTK.coeff_cost_freq([syst1, syst2, syst3], [sys1, sys2, sys3], freqs, vars)
@test cost(x0) ≈ 0 atol=sqrt(eps()) # test that cost is zero with true systems
