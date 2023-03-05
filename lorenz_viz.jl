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


#run simulation:
sim_list, markov_chain = new_run_sim(;runs=100000, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=0, rho_start=32)

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


#
new_mc_32 = []
for state in markov_chain
    if state <= 4
        push!(new_mc_32, 1)
    elseif state <= 8
        push!(new_mc_32, 2)
    else
        push!(new_mc_32, 3)
    end
end

new_ht = holding_times(new_mc, 3; dt=0.01)
new_ht_26 = holding_times(new_mc_26, 3; dt=0.01)
new_ht_32 = holding_times(new_mc_32, 3; dt=0.01)
begin
    fig = Figure(resolution=(1600, 1200))
    for i in 1:3
        ax = Axis(fig[i,1], title="$i (ρ=26)") 
        hist!(ax, new_ht_26[i])
    end
    for i in 1:3
        ax = Axis(fig[i,2], title="$i (ρ=28)") 
        hist!(ax, new_ht[i])
    end
    for i in 1:3
        ax = Axis(fig[i,3], title="$i (ρ=32)") 
        hist!(ax, new_ht_32[i])
    end
    # fig[0, 1] = Label(fig, "ρ=26")
    xlim!(fig[2,3], (0,4))
    fig
end


