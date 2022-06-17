# The ODE model.
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

params = [w1=>1, w2=>0.99, w3=>1.0, w4=>1.0, w5=>1.0, w6=>1.0, tau=>1.0, ymax=>1.0]
u0 = [A=>0.0, B=>0.0, C=>0.0, D=>0.0, E=>0.0]
tspan = (0.0, 110.0)
prob = ODEProblem(sys, u0, tspan, params)
sol = solve(prob,BS3())

plot(sol, xlabel="Time", ylabel="Activation", title="A+B Activation")



