using ModelingToolkit, ControlSystemsMTK

# Example 2 from the reference in the docstring
function testmodel(; name=:testmodel)
    @variables x1(t)=0 x2(t)=0 u(t)=0 y(t)=0
    eqs = [
        D(x1) ~ x1^2 + x2 + u
        D(x2) ~ -x2 + x1^2
        y ~ -1 - x1 + x1^2
    ]
    ODESystem(eqs, t; name)
end

model = testmodel()
@named controller = ESC(k=10, tau=0.1, a=10, w=100, wh=1000)

connections = [
    model.y ~ controller.y
    controller.u ~ model.u
]

@named closed_loop = ODESystem(connections, t, systems=[model, controller])

sys = structural_simplify(closed_loop)

x0 = [
    # model.x1 => 0.5
    # model.x2 => 0.25
    # controller.uh => -0.5
    controller.v => 0
]

prob = ODEProblem(sys, x0, (0, 10.0))

sol = solve(prob, Tsit5())

# plot(sol, vars=[model.x1, model.x2, model.y, controller.u, controller.uh], layout=5)
# display(current())

@test mean(sol[controller.uh][end-200:end]) ≈ -0.5 rtol=0.1
@test mean(sol[model.y][end-200:end]) ≈ -1.25 rtol=0.1
@test mean(sol[model.x1][end-200:end]) ≈ 0.5 rtol=0.1
@test mean(sol[model.x2][end-200:end]) ≈ 0.255 rtol=0.1
