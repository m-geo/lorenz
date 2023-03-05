using Statistics
using LinearAlgebra
# using Plots
using CairoMakie
using Latexify
using MarkovChainHammer.TransitionMatrix: perron_frobenius, generator, holding_times
using JLD
# using ColorSchemes
using HypothesisTests

include("thresholds.jl")
include("sim_utils.jl")

dt = 0.01

static_generators = Dict()
for rho in 26:1:32
    sim_list, markov_chain = new_run_sim(;runs=1e7, timing=true, 
        thresh_func=twelve_state_just_high, delta_rho=0, rho_start=rho) 
    gen = generator(markov_chain, 12; dt=0.01)
    static_generators[rho] = gen
end

# for rho in 27.5:1:30.5
#     sim_list, markov_chain = new_run_sim(;runs=1e7, timing=true, 
#         thresh_func=twelve_state_just_high, delta_rho=0, rho_start=rho) 
#     gen = generator(markov_chain, 12; dt=0.01)
#     static_generators[rho] = gen
# end

# JLD.save("static_generators_extra.jld",  "static_generators", static_generators)


sim_list_delta, markov_chain_delta = new_run_sim(;runs=10000, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=2, rho_start=26)
#
function generate_change(rho_start, delta_rho)
    changing_generators = []
    changing_mcs = []
    for i in 4:1:7
        sim_list, markov_chain = new_run_sim(;runs=10^i, timing=true, 
        thresh_func=twelve_state_just_high, delta_rho=delta_rho, rho_start=rho_start)
        gen = generator(markov_chain,12; dt=dt)
        push!(changing_generators, gen)
        push!(changing_mcs, markov_chain)
    end
    return changing_generators, changing_mcs
end

# Q_delta_bayes = BayesianGenerator(markov_chain_delta; dt = dt)
bayesian_generators = [BayesianGenerator(mc; dt=dt) for mc in changing_mcs]

function plot_rates(bayesian_generators, Q_ref; number_of_states=12)
    fig = Figure(resolution=(1600,1200))
    for i in 1:3, j in 1:4
        box = j+4*(i-1) 
        ax = Axis(fig[i,j])
        # plot!(Q_bayes.posterior.rates[box], label="1e4")
        for k in 1:3#size(bayesian_generators)
            gen = bayesian_generators[k]
            plot!(gen.posterior.rates[box], label="1e$(3+k)")
        end
        vlines!(ax, -Q_ref[box, box])
        axislegend(ax)
    end
    fig
end

plot_rates(bayesian_generators, static_generators[27])

cg_31, mcs_31 = generate_change(30, 2)

bayesians_31 = [BayesianGenerator(mc; dt=dt) for mc in mcs_31]

function plot_rates_comparison(bayesian_1, bayesian_2, Q_ref_1, Q_ref_2; number_of_states=12)
    fig = Figure(resolution=(3200,1200))
    for i in 1:3, j in 1:4
        box = j+4*(i-1) 
        ax = Axis(fig[i,j])
        for k in 1:3#size(bayesian_generators)
            gen1 = bayesian_1[k]
            plot!(gen1.posterior.rates[box], label="1e$(3+k)")
            gen2 = bayesian_2[k]
            plot!(gen2.posterior.rates[box], color="red")
        end
        vlines!(ax, [-Q_ref_1[box, box], -Q_ref_2[box, box]])
        axislegend(ax)
    end
    fig
end

plot_rates_comparison(bayesian_generators, bayesians_31, static_generators[27], static_generators[31])

function sliding_window(window_size, rho_start, rho_end; runs=1e6, dt=0.01)
    changing_generators = []
    changing_mcs = []
    # middle_values = []
    # instances = Int((rho_end-rho_start)/window_size)
    for r in rho_start:1:(rho_end-window_size)
        println("rho start is $r and rho end is $(r+window_size)")
        sim_list, markov_chain = new_run_sim(;runs=runs, timing=true, 
        thresh_func=twelve_state_just_high, delta_rho=window_size, rho_start=r)
        gen = generator(markov_chain,12; dt=dt)
        push!(changing_generators, gen)
        push!(changing_mcs, markov_chain)
        # push!(middle_values, r+window_size/2)
    end
    return changing_generators, changing_mcs#, middle_values
end

sliding_gens, sliding_mcs = sliding_window(2, 26, 32)
sliding_bayesians = [BayesianGenerator(mc; dt=dt) for mc in sliding_mcs]
sliding_reference = [static_generators[x] for x in middle_values]

sl_gen_e5, sl_mc_e5, middle_values = sliding_window(2, 26, 32; runs=1e5)
sl_bayesian_e5 = [BayesianGenerator(mc; dt=dt) for mc in sl_mc_e5]

sl_gen_e4, sl_mc_e4, middle_values = sliding_window(2, 26, 32; runs=1e4)
sl_bayesian_e4 = [BayesianGenerator(mc; dt=dt) for mc in sl_mc_e4]

middle_values = [27,28,29,30,31]

function plot_sliding_window(list_of_bayesians, ref_list; e5_list=nothing, e4_list=nothing, lims=nothing)
    fig = Figure(resolution=(3200,1200))
    colors = ["red", "orange", "green", "blue", "violet"]
    labels = [27, 28, 29, 30, 31]
    for i in 1:3, j in 1:4
        box = j+4*(i-1) 
        ax = Axis(fig[i,j])
        vlines!(ax, [-Q[box, box] for Q in ref_list])

        for k in 1:5 #generalize!!
            gen = list_of_bayesians[k]
            plot!(gen.posterior.rates[box], label="$(labels[k])", color=colors[k])
            if e5_list !== nothing
                gen = e5_list[k]
                plot!(gen.posterior.rates[box], color=colors[k], linestyle=:dash)
            end
            if e4_list !== nothing
                gen = e4_list[k]
                plot!(gen.posterior.rates[box], color=colors[k], linestyle=:dashdot)
            end
        end
        
        if lims !== nothing
            xlims!(lims[box][1],lims[box][2])
        end
        axislegend(ax)
    end
    fig
end

# lims = [(1,5),(0,30),(0,30),(1.5,5),
        # (5,7),(6,9),(6,9),(5,7),
        # (5,22),(10,80),(10,80),(7,21)]
plot_sliding_window(sliding_bayesians, sliding_reference; e5_list=sl_bayesian_e5, e4_list=sl_bayesian_e4)#, lims=lims)


# 
function kl_div(p, q)
    s1 = std(p)
    s2 = std(q)
    m1 = mean(p)
    m2 = mean(q)
    return log.(s1./s2) .+ (s1.^2 .+ (m1.-m2).^2)./(2 .* s2.^2) .- 0.5
end

kl_div(sliding_bayesians[1], sliding_bayesians[2])

function get_metrics(arg_list; test=kl_div) 
    lists = [[] for _ in 1:length(arg_list)]
    u = 0
    for bayesian_list in arg_list
        u += 1
        for dist in bayesian_list[2:end]
            # metric = kl_div(bayesian_list[1], dist)
            metric = test(dist, bayesian_list[1])
            push!(lists[u], metric)
        end
    end
    return lists
end

metrics = get_metrics((sliding_bayesians, sl_bayesian_e5, sl_bayesian_e4))#; test=HypothesisTests.ApproximateTwoSampleKSTest)

function plot_kl(metrics)
    fig = Figure(resolution=(3200,1200))
    colors = ["red", "orange", "green", "blue", "violet"]
    labels = ["1e6", "1e5", "1e4"]
    xs = [x for x in 28:31]
    for i in 1:3, j in 1:4
        box = j+4*(i-1) 
        ax = Axis(fig[i,j])
        for list_no in eachindex(metrics) #iterate 1e6 through 1e4
            #all same color + label
            tmp = Float64[]
            for k in 1:4
                push!(tmp, metrics[list_no][k][box,box])
            end
            scatter!(xs, tmp, color=colors[list_no])
            lines!(xs, tmp, color=colors[list_no], label=labels[list_no])
        end
        axislegend(ax, position=:lt)
    end
    fig
end
            

plot_kl(metrics)