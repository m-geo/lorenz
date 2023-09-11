using Statistics
using LinearAlgebra
# using Plots
using CairoMakie
using Latexify
using MarkovChainHammer.TransitionMatrix: perron_frobenius, generator, holding_times
include("thresholds.jl")
include("sim_utils.jl")

sim_list28, markov_chain28 = new_run_sim(;runs=100000, timing=true, 
        thresh_func=twelve_state_just_high, delta_rho=0, rho_start=28, init=(20,-30,100))  
#
show_full(generator(markov_chain28; dt=0.01))

generator28 = generator(markov_chain28, 12; dt=0.01)
pf28 = perron_frobenius(markov_chain28, 12)

latexify(round.(generator28, digits=3))
latexify(round.(pf28, digits=3))

## validating the smoothness of the generator changing

#run at 26, stable, for 100000 steps
#then run at 32, same. calculate Q for both, directly (generator fucn in MCH)
#then: run a long (100000) sim with full changing rho
#generate Q for subsets of that run: e.g., just the first, just the second half

sim_list26, markov_chain26 = new_run_sim(;runs=100000, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=0, rho_start=26)
sim_list32, markov_chain32 = new_run_sim(;runs=100000, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=0, rho_start=32)
sim_list_delta, markov_chain_delta = new_run_sim(;runs=100000, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=6, rho_start=26)
#
generator26 = generator(markov_chain26, 12; dt=0.01)
generator32 = generator(markov_chain32, 12; dt=0.01)
generator_delta = generator(markov_chain_delta, 12; dt=0.01)

norm(generator32-generator26) 
norm(generator_delta-generator26)
norm(generator_delta-generator32)

sim_list29, markov_chain29 = new_run_sim(;runs=100000, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=0, rho_start=29)
#
generator29 = generator(markov_chain29, 12; dt=0.01)
norm(generator29-generator26)

#halfway - calculated on first half of MC from 26 to 32, so from 26 to 29
halfway_generator = generator(markov_chain_delta[1:50000], 12; dt=0.01)
norm(generator29-generator26)
norm(halfway_generator-generator26)
norm(halfway_generator-generator29)

#so i'm concluding that it's all part of the plan. nice

#next: check the delta Q/ delta rho thing

dQ1 = generator29-generator26 # or should this be with *observed estimates* of Q?
dQ2 = generator32-generator29
norm(dQ1)
norm(dQ2)
norm(dQ2 - dQ1)
# so I guess really it's about expecting them to be the same. 
# --> I really still don't have a proper grasp on how to measure "sameness" of matrices
#but also can then vary the delta rho, see if the ratio is consisent:

norm(generator_delta - generator26)/6
norm(halfway_generator - generator26)/3

quarter_generator =  generator(markov_chain_delta[1:25000], 12; dt=0.01)
norm(quarter_generator - generator26)/1.5

three_quarter_generator = generator(markov_chain_delta[1:75000], 12; dt=0.01)
norm(three_quarter_generator - generator26)/4.5

#####  UQ for ratios #############
dt = 0.01
using MarkovChainHammer.BayesianMatrix
using Random
Random.seed!(1234)

number_of_states = 12
prior = uninformative_prior(number_of_states)
number_of_samples = 100

# mean(Q_delta_bayes_list) - mean(Q_delta_bayes)
# var(Q_delta_bayes_list) - var(Q_delta_bayes)

Q_delta_bayes = BayesianGenerator(markov_chain_delta, prior; dt = dt)
Q_delta_bayes_list = rand(Q_delta_bayes, number_of_samples)

Q_26_bayes = BayesianGenerator(markov_chain26, prior; dt=dt)
Q_26_bayes_list = rand(Q_26_bayes, number_of_samples)

full_norm_list = norm.(Q_delta_bayes_list - Q_26_bayes_list)
mean(full_norm_list./6)
std(full_norm_list./6)

Q_halfway_bayes = BayesianGenerator(markov_chain_delta[1:50000], prior; dt=dt)
Q_halfway_bayes_list = rand(Q_halfway_bayes, number_of_samples)
halfway_norm_list = norm.(Q_halfway_bayes_list - Q_26_bayes_list)
mean(halfway_norm_list./3)
std(halfway_norm_list./3)

Q_quarter_bayes = BayesianGenerator(markov_chain_delta[1:25000], prior; dt=dt)
Q_quarter_bayes_list = rand(Q_quarter_bayes, number_of_samples)
quarter_norm_list = norm.(Q_quarter_bayes_list - Q_26_bayes_list)
mean(quarter_norm_list./1.5)
std(quarter_norm_list./1.5)

## UQ for differences ####

Q_32_bayes = BayesianGenerator(markov_chain32, prior; dt=dt)
Q_32_bayes_list = rand(Q_32_bayes, number_of_samples)

mean(norm.(Q_32_bayes_list-Q_26_bayes_list))
mean(norm.(Q_delta_bayes_list-Q_26_bayes_list))
mean(norm.(Q_delta_bayes_list-Q_32_bayes_list))

Q_29_bayes = BayesianGenerator(markov_chain29, prior; dt=dt)
Q_29_bayes_list = rand(Q_29_bayes, number_of_samples)

mean(norm.(Q_29_bayes_list-Q_26_bayes_list))
mean(norm.(Q_halfway_bayes_list-Q_26_bayes_list))
mean(norm.(Q_halfway_bayes_list-Q_29_bayes_list))

###### UQ plots #####

function generate_rate_UQ(prior; number_of_samples=100, 
    markov_chain_ref=markov_chain26, markov_chain=markov_chain_delta)
    dt = 0.01
    single_ref = generator(markov_chain_ref, 12; dt=dt)
    bayesian_ref = BayesianGenerator(markov_chain_ref, prior; dt=dt)
    bayesian_ref_list = rand(bayesian_ref, number_of_samples)

    out_singles = Float64[]
    out_bayesian_means = Float64[]
    out_bayesian_stds = Float64[]

    for i in 0.1:0.1:1
        single_instance = generator(markov_chain[1:floor(Int,100000*i)], 12; dt=dt)
        bayesian_generator = BayesianGenerator(markov_chain[1:floor(Int,100000*i)], prior; dt=dt)
        bayesian_list = rand(bayesian_generator, number_of_samples)
        denom = 6*i

        push!(out_singles, norm(single_instance-single_ref)/denom)
        push!(out_bayesian_means, mean(norm.(bayesian_list-bayesian_ref_list)./denom))
        push!(out_bayesian_stds, std(norm.(bayesian_list-bayesian_ref_list)./denom))
    end
    return out_singles, out_bayesian_means, out_bayesian_stds
end

singles, bayesian_means, bayesian_sts = generate_rate_UQ(prior)
xs = Vector(0.1:0.1:1)

begin
    fig=Figure(resolution=(800,800))
    ax = Axis(fig[1,1], xlabel="progress through simulation", ylabel="norm(Q_delta-Q_26)/delta_rho")
    scatter!(xs, singles, label="single-instance", showlegend=true)
    errorbars!(xs, bayesian_means, bayesian_sts)
    scatter!(xs, bayesian_means, label="bayesian")
    ylims!(0,50)
    # Legend(fig[1,1],ax, )
    fig
end

#

function generate_difference_UQ(prior; number_of_samples=100, markov_chain_ref_start=markov_chain26,
                             markov_chain_ref_end=markov_chain32, markov_chain=markov_chain_delta)
    dt = 0.01
    single_ref_start = generator(markov_chain_ref_start, 12; dt=dt)
    single_ref_end = generator(markov_chain_ref_end, 12; dt=dt)
    bayesian_ref_start = BayesianGenerator(markov_chain_ref_start, prior; dt=dt)
    bayesian_ref_list_start = rand(bayesian_ref_start, number_of_samples)
    bayesian_ref_end = BayesianGenerator(markov_chain_ref_end, prior; dt=dt)
    bayesian_ref_list_end = rand(bayesian_ref_end, number_of_samples)

    function calculate_statistic(x, ref_start, ref_end)
        return norm(x-ref_start)+norm(x-ref_end)-norm(ref_end-ref_start)
    end

    out_singles = Float64[]
    out_bayesian_means = Float64[]
    out_bayesian_stds = Float64[]

    for i in 0.1:0.1:1
        single_instance = generator(markov_chain[1:floor(Int,100000*i)], 12; dt=dt)
        bayesian_generator = BayesianGenerator(markov_chain[1:floor(Int,100000*i)], prior; dt=dt)
        bayesian_list = rand(bayesian_generator, number_of_samples)

        push!(out_singles, calculate_statistic(single_instance, single_ref_start, single_ref_end))
        push!(out_bayesian_means, mean(calculate_statistic.(bayesian_list, bayesian_ref_list_start, bayesian_ref_list_end)))
        push!(out_bayesian_stds, std(calculate_statistic.(bayesian_list, bayesian_ref_list_start, bayesian_ref_list_end)))
    end
    return out_singles, out_bayesian_means, out_bayesian_stds
end

diff_singles, diff_means, diff_stds = generate_difference_UQ(prior)



begin
    xs = Vector(26.6:0.6:32)
    fig=Figure(resolution=(800,800))
    ax = Axis(fig[1,1], xlabel="progress through simulation", ylabel="difference")
    scatter!(xs, diff_singles, label="single-instance", showlegend=true)
    errorbars!(xs, diff_means, diff_stds)
    scatter!(xs, diff_means, label="bayesian")
    ylims!(-10,50)
    axislegend(ax)
    fig
end