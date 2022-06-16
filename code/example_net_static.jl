using DifferentialEquations
using Sundials
using ModelingToolkit
using Plots
Plots.gr(lw=2)

# fact is a normalized-Hill activation function
function act(X, W, n, EC_50)
    beta = ((EC_50^n) - 1)/(2*(EC_50^n) - 1)
    K = (beta-1)^(1/n)
    
    return W*((beta*(X^n))/((K^n)+(X^n)))
end

inhib(act_n_hill_fn, W) = W - act_n_hill_fn

AND(a, b) = a*b

OR(a, b) = a + b - a*b

@parameters max_A, tau_A, 
max_B, tau_B, 
max_C, tau_C, 
max_D, tau_D, 
max_E, tau_E, 
r1_W, r1_n, r1_EC_50, # => A
r2_W, r2_n, r2_EC_50, # => B
r3_W, r3_n, r3_EC_50, # A => C
r4_W, r4_n, r4_EC_50, # B => D
r5_W, r5_n, r5_EC_50,  # C & !D => E
r6_W, r6_n, r6_EC_50  # E => C


@variables t, A(t), B(t), C(t), D(t), E(t), 
r3_fact_A(t), r4_fact_B(t), r5_fact_C(t), r5_fact_D(t), r6_fact_E(t), r5_finhib_D(t)

Dt = Differential(t)

eqsFull = [ r3_fact_A ~ act(A, r3_W, r3_n, r3_EC_50),
        r4_fact_B ~ act(B, r4_W, r4_n, r4_EC_50),
        r5_fact_C ~ act(C, r5_W, r5_n, r5_EC_50),
        r5_fact_D ~ act(D, r5_W, r5_n, r5_EC_50),
        r5_finhib_D ~ inhib(r5_fact_D, r5_W),
        r6_fact_E ~ act(E, r6_W, r6_n, r6_EC_50),
        Dt(A) ~ (r1_W*max_A - A)/tau_A,
        Dt(B) ~ (r2_W*max_B - B)/tau_B,
        Dt(C) ~ (OR(r3_W*r3_fact_A, r6_W*r6_fact_E)*max_C - C)/tau_C,
        Dt(D) ~ (r4_W*r4_fact_B*max_D - D)/tau_D,
        Dt(E) ~ (r5_W^2*AND(r5_fact_C, r5_finhib_D)*max_E - E)/tau_E
]

@named fullsys = ODESystem(eqsFull)
fullSys = structural_simplify(fullsys)

#= 
# NOTE: mentioned in article, help for tuning
#
# From these three constraints, 
#
#   fact(0) = 0, fact(EC50) = 0.5 and fact(1) = 1.
#
# derived --> beta = (EC_50^n - 1)/(2*EC_50^n -1), K = (beta-1)^(1/n)
#
# further constrained fact(X) = 1 for X ≥ 1 to ensure that species activities are limited to YMAX.
#     As default parameters, we used W = 1, EC50 = 0.5, n = 1.4, τ = 1, and YMAX = 1.
#
# τ is the time constant for a given species
# W is the reaction weight (constrained to 0 ≤ W ≤ 1)
# YMAX is the maximal fractional activation allowing simulations of knock-down (YMAX < 1) or overexpression (YMAX > 1). (protein expression)
# Typical default reaction and species parameter values are W = 1, EC50 = 0.5, n = 1.4, τ = 1, and YMAX = 1
=#

u0 = [A => 1.0, B => 1.0, C => 0.0, D => 0.0, E => 0.0]
tend = 110.0

params = [
    max_A => 1.0, tau_A => 1.0,
    max_B => 1.0, tau_B => 1.0, 
    max_C => 1.0, tau_C => 1.0,
    max_D => 1.0, tau_D => 1.0,
    max_E => 1.0, tau_E => 1.0,
    r1_W => 0.99, r1_n => 1.4, r1_EC_50 => 0.5, # => A
    r2_W => 0.99, r2_n => 1.4, r2_EC_50 => 0.5, # => B
    r3_W => 1.0, r3_n => 1.4, r3_EC_50 => 0.5, # A => C
    r4_W => 1.0, r4_n => 1.4, r4_EC_50 => 0.5, # B => D
    r5_W => 1.0, r5_n => 1.4, r5_EC_50 => 0.5,  # C & !D => E
    r6_W => 1.0, r6_n => 1.4, r6_EC_50 => 0.5  # E => C
]


prob = ODEProblem(fullSys, u0, tend, params)

sol1 = solve(prob, QNDF(), Dabstol=1e-10, reltol=1e-10)

plot(sol1, xticks = 0:10:tend, xlims=(0.0, 50.0), ylims=(0.0, 1.1),
     title="Example_Net ( input A + B )",
     xlabel="Time (sec)",
     ylabel="Fractioanl activation")