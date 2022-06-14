using DifferentialEquations
using Sundials
using ModelingToolkit
using Plots
Plots.gr(lw=2)

# fact is a normalized-Hill activation function
function act(X, W, n, EC_50)
    beta = ((EC_50^n) - 1)/(2*(EC_50^n) -1)
    K = (beta-1)^(1/n)
    
    return W*((beta*(X^n))/((K^n)+(X^n)))
end

inhib(act_n_hill_fn, W) = W - act_n_hill_fn

AND(a, b) = a*b

OR(a, b) = a + b - a*b


# create ODE System using ModelingToolkit
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
            r6_fact_E ~ act(E, r6_W, r6_n, r6_EC_50),
            r5_finhib_D ~ inhib(r5_fact_D, r5_W),
            Dt(A) ~ (r1_W*max_A - A)/tau_A,
            Dt(B) ~ (r2_W*max_B - B)/tau_B,
            Dt(C) ~ (OR(r3_W*r3_fact_A, r6_W*r6_fact_E)*max_C - C)/tau_C,
            Dt(D) ~ (r4_W*r4_fact_B*max_D - D)/tau_D,
            Dt(E) ~ (r5_W*AND(r5_fact_C, r5_finhib_D)*max_E - E)/tau_E
]

# assign full ODESystem to fix parameter index/position
@named fullsys = ODESystem( eqsFull, t, 
    [A, B, C, D, E],
    [ max_A, tau_A,
      max_B, tau_B,
      max_C, tau_C,
      max_D, tau_D,
      max_E, tau_E,
      r1_W, r1_n, r1_EC_50,  # => A
      r2_W, r2_n, r2_EC_50,  # => B
      r3_W, r3_n, r3_EC_50,  # A => C
      r4_W, r4_n, r4_EC_50,  # B => D
      r5_W, r5_n, r5_EC_50,  # C & !D => E
      r6_W, r6_n, r6_EC_50 ] # E => C
)

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


# initial condition
u0 = [A => 0.0, B => 0.0, C => 0.0, D => 0.0, E => 0.0]

# time range
tspan = (0.0, 110.0)

# parameters
# 
# NOTE: if weight is to small, may cause solver difficult to solve answer
# by testing, choose weight = 0.99(ON_state), 0.02(OFF_state)
# 
params = [
    # species
    max_A => 1.0, tau_A => 1.0,
    max_B => 1.0, tau_B => 1.0,
    max_C => 1.0, tau_C => 1.0,
    max_D => 1.0, tau_D => 1.0,
    max_E => 1.0, tau_E => 1.0,
    # reactions
    r1_W => 0.02, r1_n => 1.4, r1_EC_50 => 0.5,  # r1 = [ => A ]
    r2_W => 0.02, r2_n => 1.4, r2_EC_50 => 0.5,  # r2 = [ => B ]
    r3_W => 1.0, r3_n => 1.4, r3_EC_50 => 0.5,   # r3 = [ A => C ]
    r4_W => 1.0, r4_n => 1.4, r4_EC_50 => 0.5,   # r4 = [ B => D ]
    r5_W => 1.0, r5_n => 1.4, r5_EC_50 => 0.5,   # r5 = [ C & !D => E ]
    r6_W => 1.0, r6_n => 1.4, r6_EC_50 => 0.5    # r6 = [ E => C ]
]

#= 
LOOKUP_TABLE：

    --- state variables u[] ---
    
    A(1), B(2), C(3), D(4), E(5)

    --- parameter p[] ---

    %%%% species %%%%
    # max_A(1), tau_A(2)
    # max_B(3), tau_B(4)
    # max_C(5), tau_C(6)
    # max_D(7), tau_D(8)
    # max_E(9), tau_E(10)

    %%%% reactions %%%%
    # r1_W(11) , r1_n(12), r1_EC_50(13)
    # r2_W(14) , r2_n(15), r2_EC_50(16)
    # r3_W(17) , r3_n(18), r3_EC_50(19)
    # r4_W(20) , r4_n(21), r4_EC_50(22)
    # r5_W(23) , r5_n(24), r5_EC_50(25)
    # r6_W(26) , r6_n(27), r6_EC_50(28)
    
    ------------------
=#

# create dynamic condition (callback method), ref = https://diffeq.sciml.ai/stable/features/callback_functions/

dosetime = [10.0, 30.0, 50.0, 70.0, 90.0] # A_on, A_off, B_on, B_off, A_&_B_on

# A_on ( t==10.0 )： A => 1.0, r1_W => 0.99
condition_1(u,t,integrator) = t == dosetime[1]
function affect_1!(integrator)
    integrator.u[1] = 1.0
    integrator.p[11] = 0.99
end
cb_1 = DiscreteCallback(condition_1, affect_1!)


# A_off ( t==30.0 )： A => 0.0, r1_W => 0.02
condition_2(u,t,integrator) = t == dosetime[2]
function affect_2!(integrator)
    integrator.u[1] = 0.0
    integrator.p[11] = 0.02
end
cb_2 = DiscreteCallback(condition_2, affect_2!)


# B_on ( t==50.0 )： B => 1.0, r2_W => 0.99
condition_3(u,t,integrator) = t == dosetime[3]
function affect_3!(integrator)
    integrator.u[2] = 1.0
    integrator.p[14] = 0.99
end
cb_3 = DiscreteCallback(condition_3, affect_3!)


# B_off ( t==70.0 )： B => 0.0, r2_W => 0.02
condition_4(u,t,integrator) = t == dosetime[4]
function affect_4!(integrator)
    integrator.u[2] = 0.0
    integrator.p[14] = 0.02
end
cb_4 = DiscreteCallback(condition_4, affect_4!)

# A_B_on ( t==90.0 )： A => 1.0, r1_W => 0.99, B => 1.0, r2_W => 0.99
condition_5(u,t,integrator) = t == dosetime[5]
function affect_5!(integrator)
    integrator.u[1] = 1.0
    integrator.p[11] = 0.99
    integrator.u[2] = 1.0
    integrator.p[14] = 0.99
end
cb_5 = DiscreteCallback(condition_5, affect_5!)

# collect all callback functions
cb_set = CallbackSet(cb_1, cb_2, cb_3, cb_4, cb_5)

# initial ODE Problem
prob = ODEProblem(fullSys, u0, tspan, params)

# solve problem
@time sol1 = solve(prob, QNDF(), callback=cb_set, tstops=dosetime)


plot(sol1, xticks = 0:5:110, yticks = 0:0.1:1, xlims=(0, 130), ylims=(0.0, 1.1),
     title="Example_Net (Dynamic reaction)", legend=:bottomright, 
     xlabel="Time (sec)",
     ylabel="Fractioanl activation")