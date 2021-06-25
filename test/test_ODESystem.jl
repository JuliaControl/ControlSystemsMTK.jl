using ControlSystemsMTK, ControlSystems, ModelingToolkit, DifferentialEquations, RobustAndOptimalControl
## Test SISO (single input, single output) system
@parameters t

P0 = tf(1.0, [1, 1])                         |> ss
C0 = pid(kp = 1, ki = 1) * tf(1, [0.01, 1])  |> ss

@named P           = ODESystem(P0)
# @named nonlinear_P = sconnect(x->sign(x)*sqrt(abs(x)), P) # apply input-nonlinearity
@named C           = ODESystem(C0)
@named loopgain    = sconnect(C, P)
@named fb          = feedback(loopgain, sin(t))
@show fb           = structural_simplify(fb)

@test length(states(P)) == 3 # 1 + u + y
@test length(states(C)) == 4 # 2 + u + y

x0 = Pair[
    # fb.loopgain.nonlinear_P.P.x1 => 1 # a bit inconvenient to specify initial states
    # loopgain.nonlinear_P.P.x1 => 1
    loopgain.P.x1 => 1
    # P.x1 => 1
]
p = Pair[]

prob = ODEProblem(fb, x0, (0.0, 10.0), p)
sol = solve(prob, Tsit5())   
plot(sol)
# xtraj = Array(sol)[3,:]
# lsim(fb, )


# === Go the other way, ODESystem -> StateSpace ================================
x = states(P) # I haven't figured out a good way to access states, so this is a bit manual and ugly
@nonamespace P02_named = named_ss(P, P.u, [P.y])
@test P02_named.x_names == [Symbol("x1(t)")]
@test P02_named.u_names == [Symbol("u(t)")]
@test P02_named.y_names == [Symbol("y(t)")]

@nonamespace P02 = ss(P, P.u, P.y)
@test P02 == P0 # verify that we get back what we started with

# same for controller
x = states(C) 
@nonamespace C02 = ss(C, [C.u], C.y)
@test C02 == C0

# Now try do the same with the feedback interconnection. This fails, I cannot figure out how to provide the input I want to linearize w.r.t. to. Below is an attempt by creating an input variable `r`, but that causes `structural_simplify` to complain.
@register r(t)
@register rfun(t)
@named fbu = feedback(loopgain, rfun(t))
@show fbu = structural_simplify(fbu)
fbus = states(fbu)
fb2 = @nonamespace ss(fbu, rfun(t), fbu.y)
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

P0 = ss(A,B,C,0)

s = tf("s")
Cv0 = tx * pid(kp = 4, ki = 2) * tf(1, [0.01, 1])  #|> ss
Cp0 = tx * pid(kp = 2)         * tf(1, [0.01, 1])  #|> ss

P0i = deepcopy(P0)
P0i.D .= 1e-8

@named RF = ODESystem(tf(1, [0.001, 1])) # ref filter
@named RFv = ODESystem(tf(1, [0.001, 1]))
@named P  = ODESystem(P0) # system model
@named Fv = ODESystem(inv(P0i[2,1])) # feedforward vel
@named Fp = ODESystem(inv(P0i[1,1])) # feedforward pos
@named Cv = ODESystem(Cv0) # vel controller
@named Cp = ODESystem(Cp0) # pos controller
fb = let ref0 = 0.2sin(t), disturbance = 100sign(t-10)
    @variables u(t) y(t) 
    Dₜ = Differential(t)
    fb = ODESystem([
        u    ~ ref0 # The input 
        ref0 ~ RF.u # filter position reference
        expand_derivatives(Dₜ(ref0) ~ RFv.u) # Filter vel reference
        RF.y - P.y1 ~ Cp.u # pos controller input is pos error
        disturbance ~ P.u2 # disturbance enters on arm acceleration
        Cp.y + RFv.y - P.y2 ~ Cv.u # vel controller input is vel error + pos controller output
        Cv.y + Fp.y ~ P.u1 # robot input is vel controller output and trq ff
        y ~ P.y1 # output is robot motor pos
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
plot!(sol.t, sol[P.x1]-sol[P.x3], sp=12, lab="Δq")

##
plot(sol.t, sol[P.x1]-@nonamespace(sol[fb.u]), lab="qₘ", title="Control error")
plot!(sol.t, sol[P.x3]-@nonamespace(sol[fb.u]), lab="qₐ")