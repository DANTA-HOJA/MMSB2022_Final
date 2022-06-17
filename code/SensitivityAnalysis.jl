using OrdinaryDiffEq, CairoMakie
using ModelingToolkit
using LinearAlgebra
using DifferentialEquations
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


# create ODE System using ModelingToolkit

# --- species ---
@parameters( 
    max_A,
    max_B, 
    max_C, 
    max_D,
    max_E,
    tau,
    n,
    EC_50)
# --- reactions ---
@parameters(
    r1_W,   # r1 = [ => A ]
    r2_W,  # r2 = [ => B ]
    r3_W,  # r3 = [ A => C ]
    r4_W,  # r4 = [ B => D ]
    r5_W,   # r5 = [ C & !D => E ]
    r6_W )  # r6 = [ E => C ]  

@variables t, A(t), B(t), C(t), D(t), E(t), 
r3_fact_A(t), r4_fact_B(t), r5_fact_C(t), r5_fact_D(t), r6_fact_E(t), r5_finhib_D(t)

Dt = Differential(t)

eqsFull = [ r3_fact_A ~ act(A, r3_W, n, EC_50),
            r4_fact_B ~ act(B, r4_W, n, EC_50),
            r5_fact_C ~ act(C, r5_W, n, EC_50),
            r5_fact_D ~ act(D, r5_W, n, EC_50),
            r6_fact_E ~ act(E, r6_W, n, EC_50),
            r5_finhib_D ~ inhib(r5_fact_D, r5_W),
            Dt(A) ~ (r1_W*max_A - A)/tau,
            Dt(B) ~ (r2_W*max_B - B)/tau,
            Dt(C) ~ (OR(r3_W*r3_fact_A, r6_W*r6_fact_E)*max_C - C)/tau,
            Dt(D) ~ (r4_W*r4_fact_B*max_D - D)/tau,
            Dt(E) ~ (r5_W*AND(r5_fact_C, r5_finhib_D)*max_E - E)/tau
]

# assign full ODESystem to fix parameter index/position
@named fullsys = ODESystem(eqsFull, t, 
    [A, B, C, D, E],
    [ max_A, 
      max_B, 
      max_C, 
      max_D, 
      max_E, 
      tau,
      n,
      EC_50,
      r1_W,  # r1 = [ => A ]
      r2_W,  # r2 = [ => B ]
      r3_W, # r3 = [ A => C ]
      r4_W,  # r4 = [ B => D ]
      r5_W,  # r5 = [ C & !D => E ]
      r6_W ] # r6 = [ E => C ]
)

sys = structural_simplify(fullsys)


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


# initial condition
u0 = [A => 0.0, B => 0.0, C => 0.0, D => 0.0, E => 0.0]

# time range
tspan = (0.0, 100.0)

# parameters
# 
# NOTE: if weight is to small, may cause solver difficult to solve answer
# by testing, choose weight = 0.99(ON_state), 0.02(OFF_state)
# 
params = [
    # species
    max_A => 1.0, 
    max_B => 1.0, 
    max_C => 1.0, 
    max_D => 1.0, 
    max_E => 1.0, 
    tau   => 1.0,
    n     => 1.4,
    EC_50 => 0.5,
    # reactions
    r1_W => 0.01,  # r1 = [ => A ]
    r2_W => 0.01,  # r2 = [ => B ]
    r3_W => 1.0,   # r3 = [ A => C ]
    r4_W => 1.0,   # r4 = [ B => D ]
    r5_W => 1.0,   # r5 = [ C & !D => E ]
    r6_W => 1.0    # r6 = [ E => C ]
]

params1 = [
    # species
    max_A => 0.75, 
    max_B => 1.0, 
    max_C => 1.0, 
    max_D => 1.0, 
    max_E => 1.0, 
    tau   => 1.0,
    n     => 1.4,
    EC_50 => 0.5,
    # reactions
    r1_W => 0.01,  # r1 = [ => A ]
    r2_W => 0.01,  # r2 = [ => B ]
    r3_W => 1.0,   # r3 = [ A => C ]
    r4_W => 1.0,   # r4 = [ B => D ]
    r5_W => 1.0,   # r5 = [ C & !D => E ]
    r6_W => 1.0    # r6 = [ E => C ]
]
params2 = [
    # species
    max_A => 1.0, 
    max_B => 0.75, 
    max_C => 1.0, 
    max_D => 1.0, 
    max_E => 1.0, 
    tau   => 1.0,
    n     => 1.4,
    EC_50 => 0.5,
    # reactions
    r1_W => 0.01,  # r1 = [ => A ]
    r2_W => 0.01,  # r2 = [ => B ]
    r3_W => 1.0,   # r3 = [ A => C ]
    r4_W => 1.0,   # r4 = [ B => D ]
    r5_W => 1.0,   # r5 = [ C & !D => E ]
    r6_W => 1.0    # r6 = [ E => C ]
]
params3 = [
    # species
    max_A => 1.0, 
    max_B => 1.0, 
    max_C => 0.75, 
    max_D => 1.0, 
    max_E => 1.0, 
    tau   => 1.0,
    n     => 1.4,
    EC_50 => 0.5,
    # reactions
    r1_W => 0.01,  # r1 = [ => A ]
    r2_W => 0.01,  # r2 = [ => B ]
    r3_W => 1.0,   # r3 = [ A => C ]
    r4_W => 1.0,   # r4 = [ B => D ]
    r5_W => 1.0,   # r5 = [ C & !D => E ]
    r6_W => 1.0    # r6 = [ E => C ]
]
params4 = [
    # species
    max_A => 1.0, 
    max_B => 1.0, 
    max_C => 1.0, 
    max_D => 0.75, 
    max_E => 1.0, 
    tau   => 1.0,
    n     => 1.4,
    EC_50 => 0.5,
    # reactions
    r1_W => 0.01,  # r1 = [ => A ]
    r2_W => 0.01,  # r2 = [ => B ]
    r3_W => 1.0,   # r3 = [ A => C ]
    r4_W => 1.0,   # r4 = [ B => D ]
    r5_W => 1.0,   # r5 = [ C & !D => E ]
    r6_W => 1.0    # r6 = [ E => C ]
]
params5 = [
    # species
    max_A => 1.0, 
    max_B => 1.0, 
    max_C => 1.0, 
    max_D => 1.0, 
    max_E => 0.75, 
    tau   => 1.0,
    n     => 1.4,
    EC_50 => 0.5,
    # reactions
    r1_W => 0.01,  # r1 = [ => A ]
    r2_W => 0.01,  # r2 = [ => B ]
    r3_W => 1.0,   # r3 = [ A => C ]
    r4_W => 1.0,   # r4 = [ B => D ]
    r5_W => 1.0,   # r5 = [ C & !D => E ]
    r6_W => 1.0    # r6 = [ E => C ]
]

prob = ODEProblem(sys, u0, tspan, params)
sol = solve(prob,QNDF())
y_orig = sol.u[length(sol.u),:][1]
#y_orig = mapreduce(permutedims, hcat, transpose.(y_fin_orig))

prob1= ODEProblem(sys, u0, tspan, params1)
prob2= ODEProblem(sys, u0, tspan, params2)
prob3= ODEProblem(sys, u0, tspan, params3)
prob4= ODEProblem(sys, u0, tspan, params4)
prob5= ODEProblem(sys, u0, tspan, params5)
sol1 = solve(prob1,QNDF())
sol2 = solve(prob2,QNDF())
sol3 = solve(prob3,QNDF())
sol4 = solve(prob4,QNDF())
sol5 = solve(prob5,QNDF())

y_final1 = sol1.u[length(sol1.u),:][1]
y_final2 = sol2.u[length(sol2.u),:][1]
y_final3 = sol3.u[length(sol3.u),:][1]
y_final4 = sol4.u[length(sol4.u),:][1]
y_final5 = sol5.u[length(sol5.u),:][1]

#y_fin = mapreduce(permutedims, hcat, transpose.(y_fin)) 
#y_final =vcat(y_final, y_fin)      # Make column of y_final equal to steady state values of y_new

# Generate normalized sensitivity coefficients: S=dY/dp*p/Y
function sensitivity(y_final, y_orig, dP, p0)
    dY= y_final - y_orig
    S = (dY ./ dP) .* p0 ./ y_orig
    #S1 = getindex.(S,1)
    return S
end

SA = sensitivity(y_final1, y_orig, -0.25 , 1)[:,:]
SB = sensitivity(y_final2, y_orig, -0.25 , 1)[:,:]
SC = sensitivity(y_final3, y_orig, -0.25 , 1)[:,:]
SD = sensitivity(y_final4, y_orig, -0.25 , 1)[:,:]
SE = sensitivity(y_final5, y_orig, -0.25 , 1)[:,:]

S2d =hcat(SA, SB, SC, SD, SE) 
S2d

fig = Figure(resolution = (600, 400))
ax, hm = CairoMakie.heatmap(fig[1,2], S2d,
colormap = :diverging_gkr_60_10_c40_n256, colorrange=(-1.5, 1.5),
axis = (; xlabel="Pertubed species", ylabel="Sensitivity of output species"))
Colorbar(fig[1, 1], hm, ticks = -1.5:0.5:1.5)
fig


#= ensemble analysis
a = map(LinRange(0.0, 1000.0, 50)) do k1
    p = copy(params)
    p[k_1] = k1
	prob = SteadyStateProblem(sys, u0, p)
	sol = solve(prob, DynamicSS(Rodas5()))
	sol[1]
end

function prob_func(prob,i,repeat)
    x = 0.3rand(2)
    remake(prob,p=[p[1:2];x])
  end

ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
sim = solve(ensemble_prob,SRIW1(),trajectories=10)
=#