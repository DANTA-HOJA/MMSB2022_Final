# The ODE model.
using DiffEqCallbacks
using DifferentialEquations
using ModelingToolkit
using Plots
using LinearAlgebra

# plotting default option(s)
Plots.gr(linewidth=2)

# OR logic gate
function OR(x,y)
    z = x + y - x*y
    return z
end

# AND logic gate, multiplying all of the reactants together
AND(x,y) = x*y

# Convenience functions
hill(x, k) = x / (x + k)
hill(x, k, n) = hill(x^n, k^n)

# hill activation function with parameters w (weight), n (Hill coeff), K
function act(x,w)
    n = 1.4
    EC_50 = 0.5
    β = ((EC_50^n) - 1)/(2*(EC_50^n) - 1)
    K = (β - 1)^(1/n)
    fact = w * β * hill(x, K, n)
    return fact
end

# inverse hill function with parameters w (weight), n (Hill coeff), K 
function inhib(x,w)
    finhib = w - act(x,w)
    return finhib
end


#= reaction crosstalk
=> A            w1
=> B            w2
A => C          w3
B => D          w4
C & !D => E     w5
E => C          w6
=#
@parameters w1 w2 w3 w4 w5 w6 tau ymax
@variables t, A(t), B(t), C(t), D(t), E(t),
r3_fact_A(t), r4_fact_B(t), r5_fact_C(t), r5_fact_D(t), r6_fact_E(t), r5_finhib_D(t)

Diff = Differential(t)

eqs = [  r3_fact_A ~ act(A, w3),
         r4_fact_B ~ act(B, w4),
         r5_fact_C ~ act(C, w5),
         r5_fact_D ~ act(D, w5),
         r6_fact_E ~ act(E, w6),
         r5_finhib_D ~ inhib(D,  w5),
        Diff(A) ~  (w1 * ymax - A) / tau,
        Diff(B) ~  (w2 * ymax - B) / tau,
        Diff(C) ~  (OR(w3 *r3_fact_A, w6 * r6_fact_E) * ymax  - C) / tau,
        Diff(D) ~  (w4 * r4_fact_B * ymax - D) / tau,
        Diff(E) ~  (w5 *AND(r5_fact_C, r5_finhib_D) * ymax - E) / tau  
]

# assign full ODESystem to fix parameter index/position
@named fullsys = ODESystem(eqs, t, [A, B, C, D, E], [w1, w2, w3, w4, w5, w6, tau, ymax])
sys = structural_simplify(fullsys)


params = [w1=>0.0, w2=>0.0, w3=>1.0, w4=>1.0, w5=>1.0, w6=>1.0, tau=>1.0, ymax=>1.0]
u0 = [A=>0.0, B=>0.0, C=>0.0, D=>0.0, E=>0.0]
stop=[10, 30, 50, 70, 90]
tspan = (0.0, 110.0)
prob = ODEProblem(sys, u0, tspan, params)

condition_1(u,t,integrator) = t == stop[1]
function affect_1!(integrator) 
    integrator.p[1] = 0.99
    integrator.u[1] = 1.0
end
cb1 = DiscreteCallback(condition_1, affect_1!)

condition_2(u,t,integrator) = t == stop[2]
function affect_2!(integrator)
    integrator.p[1] = 0.01
    integrator.u[1] = 0
end
cb2 = DiscreteCallback(condition_2, affect_2!)

condition_3(u,t,integrator) = t == stop[3]
function affect_3!(integrator) 
    integrator.p[2] = 0.99
    integrator.u[2] = 1.0
end
cb3 = DiscreteCallback(condition_3, affect_3!)

condition_4(u,t,integrator) = t == stop[4]
function affect_4!(integrator) 
    integrator.p[2] = 0.01
    integrator.u[2] = 0
end
cb4 = DiscreteCallback(condition_4, affect_4!)

condition_5(u,t,integrator) = t ==stop[5]
function affect_5!(integrator) 
    integrator.p[2] = 0.99
    integrator.p[1] = 0.99
    integrator.u[2] = 1.0
    integrator.u[1] = 1.0
end
cb5 = DiscreteCallback(condition_5, affect_5!)

cb_set = CallbackSet(cb1, cb2, cb3, cb4, cb5)
sol = solve(prob, BS3(), tstops = stop , callback=cb_set)

plot(sol, xlabel="Time", ylabel="Activation", title="Example_Net", legend=:outerbottomright) 

# plot figure with split AB and CDE
Fig_AB = plot()
plot!(Fig_AB, sol, vars=([1,2]))

Fig_CDE = plot()
plot!(Fig_CDE, sol, vars=([3,4,5]), xlabel="Time (sec)")


