using Statistics
using LinearAlgebra
# using Plots
using CairoMakie
# using Latexify
using MarkovChainHammer.TransitionMatrix: perron_frobenius, generator, holding_times
using MarkovChainHammer.BayesianMatrix
using JLD
using Distributions
# using ColorSchemes
# using HypothesisTests

include("thresholds.jl")
include("sim_utils.jl")

dt = 0.01

# static_generators = Dict()
# for rho in 26:1:32
#     sim_list, markov_chain = new_run_sim(;runs=1e7, timing=true, 
#         thresh_func=twelve_state_just_high, delta_rho=0, rho_start=rho) 
#     gen = generator(markov_chain, 12; dt=0.01)
#     static_generators[rho] = gen
# end

# for rho in 27.5:1:30.5
#     sim_list, markov_chain = new_run_sim(;runs=1e7, timing=true, 
#         thresh_func=twelve_state_just_high, delta_rho=0, rho_start=rho) 
#     gen = generator(markov_chain, 12; dt=0.01)
#     static_generators[rho] = gen
# end

# further_generators = Dict()
# for rho in 33:1:36
#     sim_list, markov_chain = new_run_sim(;runs=1e7, timing=true, 
#         thresh_func=twelve_state_just_high, delta_rho=0, rho_start=rho) 
#     gen = generator(markov_chain, 12; dt=0.01)
#     further_generators[rho] = gen
# end

# JLD.save("./vars/further_generators.jld",  "further_generators", further_generators)
# JLD.save("static_generators_extra.jld",  "static_generators", static_generators)
static_generators = JLD.load("./vars/static_generators.jld", "static_generators")


##########

# bayesian_static_generators = Dict()
# for rho in 26:1:36
#     println("starting rho=$rho")
#     sim_list, markov_chain = new_run_sim(;runs=1e7, timing=true, 
#         thresh_func=twelve_state_just_high, delta_rho=0, rho_start=rho) 
#     gen = BayesianGenerator(markov_chain; dt=0.01)
#     bayesian_static_generators[rho] = gen
# end


# static_means = Dict()
# static_stds = Dict()
# for rho in 26:1:36
#     Q_bayes = bayesian_static_generators[rho]
#     static_means[rho] = mean(Q_bayes)
#     static_stds[rho] = std(Q_bayes)
# end
# JLD.save("./vars/bayesian_static_means.jld",  "bayesian_static_means", static_means)
# JLD.save("./vars/bayesian_static_stds.jld",  "bayesian_static_stds", static_stds)

# bayesian_static_generators = JLD.load("./vars/bayesian_static_generators.jld", "bayesian_static_generators")


sim_list_delta, markov_chain_delta = new_run_sim(;runs=100000, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=0, rho_start=26)

perron_frobenius(markov_chain_delta; step=10)

###
function generate_change(rho_start, delta_rho)
    # 
    changing_generators = []
    changing_mcs = []
    for i in 4:1:7 #iterate through degree of runs
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


### this is where we actually start?

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

middle_values = [27,28,29,30,31]

sliding_gens, sliding_mcs = sliding_window(2, 26, 32)
sliding_bayesians = [BayesianGenerator(mc; dt=dt) for mc in sliding_mcs]
sliding_reference = [static_generators[x] for x in middle_values]

#slidigng bayesians is a list of bayesian generators for each of the windows (5 total)

sl_gen_e5, sl_mc_e5 = sliding_window(2, 26, 32; runs=1e5)
sl_bayesian_e5 = [BayesianGenerator(mc; dt=dt) for mc in sl_mc_e5]

sl_gen_e4, sl_mc_e4 = sliding_window(2, 26, 32; runs=1e4)
sl_bayesian_e4 = [BayesianGenerator(mc; dt=dt) for mc in sl_mc_e4]


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
function kl_div(p, q; significance=nothing)
    # for significance: integer indicating factor of sigma difference
    s1 = std(p)
    s2 = std(q)
    m1 = mean(p)
    m2 = mean(q)
    metric = log.(s1./s2) .+ (s1.^2 .+ (m1.-m2).^2)./(2 .* s2.^2) .- 0.5
    if significance === nothing
        return metric
    else
        thresh =  log.(s1./s2) .+ (s1.^2 .+ (significance .* s1).^2)./(2 .* s2.^2) .- 0.5
        return metric .>= thresh
    end
end

kl_div(sliding_bayesians[1], sliding_bayesians[2])

function get_metrics(arg_list; significance=nothing) 
    lists = [[] for _ in 1:length(arg_list)]
    signif = [[] for _ in 1:length(arg_list)]
    u = 0
    for bayesian_list in arg_list
        u += 1
        for dist in bayesian_list[2:end]
            metric = kl_div(bayesian_list[1], dist) #first the reference, then one under question
            push!(lists[u], metric)
            if significance !== nothing
                sig = kl_div(bayesian_list[1], dist; significance=significance)
                push!(signif[u], sig)
            end
        end
    end
    return lists, signif
end

metrics, metrics_sig = get_metrics((sliding_bayesians, sl_bayesian_e5, sl_bayesian_e4); significance=2)

function plot_kl(metrics; metrics_sig)
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
            tmp_sig = Float64[]
            marker_shapes = []
            for k in 1:4
                push!(tmp, metrics[list_no][k][box,box])
                # push!(tmp, metrics_sig[list_no][k][box,box]) #condition this!!!
                if metrics_sig[list_no][k][box,box] 
                    push!(marker_shapes, :circle) #o
                else
                    push!(marker_shapes, :diamond) #s
                end
            end
            # marker_shapes = [tmp_sig[i] ? 'o' : 's' for i in 1:length(tmp_sig)]
            # print(marker_shapes)
            scatter!(xs, tmp, color=(colors[list_no], 0.5), marker=marker_shapes, markersize=20)
            lines!(xs, tmp, color=colors[list_no], label=labels[list_no])
        end
        axislegend(ax, position=:lt)
    end
    fig
end
            
plot_kl(metrics; metrics_sig)


####### off-diagonal entries

function plot_sliding_off_diag(list_of_bayesians, ref_list; e5_list=nothing, e4_list=nothing, lims=nothing)
    fig = Figure(resolution=(800,400))
    colors = ["red", "orange", "green", "blue", "violet"]
    labels = [27, 28, 29, 30, 31]
    box_list = [(5,9), (8,12)]
    for n in eachindex(box_list)
        i, j = box_list[n]
        ax = Axis(fig[1,n], title="$i -> $j")
        # vlines!(ax, [Q[i, j] for Q in ref_list])

        for k in 1:5 #generalize!?
            gen = list_of_bayesians[k]
            alpha_j = gen.posterior.exit_probabilities[i].alpha[j-1] #j-1 only works if looking at transitions into later states!
            alpha_0 = sum(gen.posterior.exit_probabilities[i].alpha)
            dist = Distributions.Beta(alpha_j, alpha_0 - alpha_j)
            plot!(dist, label="$(labels[k])", color=colors[k])
            if e5_list !== nothing
                gen = e5_list[k]
                alpha_j = gen.posterior.exit_probabilities[i].alpha[j-1] #j-1 only works if looking at transitions into later states!
                alpha_0 = sum(gen.posterior.exit_probabilities[i].alpha)
                dist = Distributions.Beta(alpha_j, alpha_0 - alpha_j)
                plot!(dist,  color=colors[k], linestyle=:dash)
            end
            if e4_list !== nothing
                gen = e4_list[k]
                alpha_j = gen.posterior.exit_probabilities[i].alpha[j-1] #j-1 only works if looking at transitions into later states!
                alpha_0 = sum(gen.posterior.exit_probabilities[i].alpha) #double check that ? ^ like. i vs j.
                dist = Distributions.Beta(alpha_j, alpha_0 - alpha_j)
                plot!(dist,  color=colors[k], linestyle=:dashdot)
            end
        end
        
        if lims !== nothing
            xlims!(lims[n][1],lims[n][2])
        end
        axislegend(ax)
    end
    fig
end

# lims = [(0,0.75),(0,30),(0,30)]
plot_sliding_off_diag(sliding_bayesians, sliding_reference; e5_list=sl_bayesian_e5, e4_list=sl_bayesian_e4)#, lims=lims)
