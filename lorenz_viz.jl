using Statistics
# using Clustering
# using Polynomials
using LinearAlgebra
# using Makie
using Plots
using CairoMakie
# using Latexify
using MarkovChainHammer.TransitionMatrix: perron_frobenius, generator, holding_times
include("thresholds.jl")
include("sim_utils.jl")

#focus is mostly on plotting historgrams


#run simulation:
sim_list_26, markov_chain_26 = new_run_sim(;runs=1e7, timing=true, 
                            thresh_func=z_tercile_thresh, delta_rho=0, rho_start=26)
sim_list_delta, markov_chain_delta = new_run_sim(;runs=1e7, timing=true, 
                            thresh_func=z_tercile_thresh, delta_rho=6, rho_start=26)                            
sim_list_32, markov_chain_32 = new_run_sim(;runs=1e7, timing=true, 
                            thresh_func=z_tercile_thresh, delta_rho=0, rho_start=32)

## visualize holding times
ht = holding_times(markov_chain, 12; dt=0.01)

hist(ht[12])

begin
    fig = Figure(resolution=(1600, 1200))
    for i in 1:3, j in 1:4
        box = j+4*(i-1) #where 3 is from the max j
        ax = Axis(fig[i,j], title="$box") 
        hist!(ax, ht[box])
    end
    fig
end


# aggregate the states
# new_mc_32 = []
# for state in markov_chain
#     if state <= 4
#         push!(new_mc_32, 1)
#     elseif state <= 8
#         push!(new_mc_32, 2)
#     else
#         push!(new_mc_32, 3)
#     end
# end

ht_delta = holding_times(markov_chain_delta, 3; dt=0.01)
ht_26 = holding_times(markov_chain_26, 3; dt=0.01)
ht_32 = holding_times(markov_chain_32, 3; dt=0.01)
begin
    fig = Figure(resolution=(1600, 1200))
    for i in 1:3
        ax = Axis(fig[i,1], title="$i (ρ=26)") 
        hist!(ax, ht_26[i])
    end
    for i in 1:3
        ax = Axis(fig[i,2], title="$i (ρ changing)") 
        hist!(ax, ht_delta[i])
    end
    for i in 1:3
        ax = Axis(fig[i,3], title="$i (ρ=32)") 
        hist!(ax, ht_32[i])
    end
    # fig[0, 1] = Label(fig, "ρ=26")
    # xlim!(fig[2,3], (0,4))
    fig
end


