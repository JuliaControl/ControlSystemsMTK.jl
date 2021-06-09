## Test SISO (single input, single output) system
@parameters t

P0 = tf(1.0, [1, 1])                         |> ss
C0 = pid(kp = 1, ki = 1) * tf(1, [0.01, 1])  |> ss

@named P           = ODESystem(P0)
# @named nonlinear_P = connect(x->sign(x)*sqrt(abs(x)), P) # apply input-nonlinearity
@named C           = ODESystem(C0)
@named loopgain    = connect(C, P)
@named fb          = feedback(loopgain, sin(t))
@show fb = structural_simplify(fb)

x0 = Pair[
    # fb.loopgain.nonlinear_P.P.x1 => 1 # a bit inconvenient to specify initial states
    # loopgain.nonlinear_P.P.x1 => 1
    loopgain.P.x1 => 1
    # P.x1 => 1
]
p = Pair[]

prob = ODEProblem(simplified_sys, x0, (0.0, 10.0), p)
sol = solve(prob, Tsit5())
plot(sol)

# === Go the other way, ODESystem -> StateSpace ================================
using Test
x = ModelingToolkit.get_states(P) # I haven't figured out a good way to access states, so this is a bit manual and ugly
P02 = ss(P, x[2:2], x[3:3])
@test P02 == P0 # verify that we get back what we started with

# same for controller
x = ModelingToolkit.get_states(C) 
C02 = ss(C, x[3:3], x[4:4])
@test C02 == C0

# Now try do the same with the feedback interconnection. This fails, I cannot figure out how to provide the input I want to linearize w.r.t. to. Below is an attempt by creating an input variable `r`, but that causes `structural_simplify` to complain.
@variables r(t)
@named fbu = feedback(loopgain, r)
@show fbu = structural_simplify(fbu)
fbus = ModelingToolkit.get_states(fbu)
fb2 = @nonamespace ss(fbu, [r], [fbu.y])
feedback(P0*C0) # fb2 should be similar to this feeback interconnection calculated by ControlSystems




## Test DMM with cascade =====================================================
@parameters t

Jm = 10
Ja = 30
k  = 4e4
c0 = 100
c1 = 100
c2 = 0.001 # This should be small to see large oscillations on arm side. Doesn't make much sense to have damping on arm-side anyways
tx = 20

A = [
    0.0 1 0 0
    -k/Jm -(c1 + c0)/Jm k/Jm c1/Jm
    0 0 0 1
    k/Ja c1/Ja -k/Ja -(c1 + c2)/Ja
]

B  = [0, 1/Jm, 0, 0]
Bd = [0, 0, 0, 1/Ja]
B = [B Bd]
C = [
    1 0 0 0
    0 1 0 0
]

@info "Resonance frequency: $((imag.(eigvals(A)) ./ (2π)))"

#
using NNlib: σ
P0 = ss(A,B,C,0)

s = tf("s")
# estun uses parameters in the range Kp ≈ 50, Kv ≈ 300, Ki ≈ 300 (time)
Cv0 = tx * pid(kp = 4, ki = 2) * tf(1, [0.01, 1])  #|> ss
Cp0 = tx * pid(kp = 2)         * tf(1, [0.01, 1])  #|> ss

P0i = deepcopy(P0)
P0i.D .= 1e-8

@named RF = ODESystem(tf(1, [0.001, 1]))
@named RFv = ODESystem(tf(1, [0.001, 1]))
@named P  = ODESystem(P0)
@named Fv = ODESystem(inv(P0i[2,1]))
@named Fp = ODESystem(inv(P0i[1,1]))
@named Cv = ODESystem(Cv0)
@named Cp = ODESystem(Cp0)
fb = let ref0 = 0.2sin(t), disturbance = 100sign(t-10)
    @variables u(t) y(t) 
    Dₜ = Differential(t)
    fb = ODESystem([
        u    ~ ref0 # The input 
        ref0 ~ RF.u # filter position reference
        expand_derivatives(Dₜ(ref0) ~ RFv.u) # Filter vel reference
        RF.y - P.y₁ ~ Cp.u # pos controller input is pos error
        disturbance ~ P.u₂ # disturbance enters on arm acceleration
        Cp.y + RFv.y - P.y₂ ~ Cv.u # vel controller input is vel error + pos controller output
        Cv.y + Fp.y ~ P.u₁ # robot input is vel controller output and trq ff
        y ~ P.y₁ # output is robot motor pos
        RF.y ~ Fp.u
    ], t; systems=[P, Cv, Cp, Fp, RF, RFv], name=:feedback)
end
simplified_sys = structural_simplify(fb)


x0 = Pair[
    P.x1 => 0.0
    P.x3 => 0.0
]
p = Pair[]

prob = ODEProblem(simplified_sys, x0, (0.0, 20.0), p, jac=true)
sol = solve(prob, OrdinaryDiffEq.Rodas5(), rtol=1e-8, atol=1e-8, saveat=0:0.01:20)
@show sol.retcode
plot(sol, layout=length(states(simplified_sys))+1)
plot!(sol.t, sol[P.x₁]-sol[P.x₃], sp=12, lab="Δq")

##
plot(sol.t, sol[P.x₁]-@nonamespace(sol[fb.u]), lab="qₘ", title="Control error")
plot!(sol.t, sol[P.x₃]-@nonamespace(sol[fb.u]), lab="qₐ")