using Statistics
# using Clustering
# using Polynomials
using LinearAlgebra
# using Makie
# using StatsBase
using Plots
using Latexify
include("thresholds.jl")
include("sim_utils.jl")

function show_full(array)
    show(IOContext(stdout, :limit=>false), MIME"text/plain"(), array)
end

#run simulation:
X, Y, Z, state_list, holding_times = run_sim(;runs=100000, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=0, rho_start=28)
number_of_states = 12

##
Plots.histogram(Z, legend=:none, normalize=true)
Plots.vline!([22.5, 37.90], linewidth=4, legend=:none)
Plots.xlabel!("Z")
Plots.ylabel!("frequency")
Plots.savefig("z_histogram.png")

sim_list, mc_chain = new_run_sim(;runs=100000, timing=true, 
    thresh_func=twelve_state_just_high, delta_rho=3, rho_start=26)

points = transpose(hcat(X, Y, Z))
# points_list = [line for row in [X Y Z][1:300,:]]
##

## construct probability matrix directly
norm_constr_prob, empirical_T = construct_p(number_of_states, state_list, holding_times, simple=false)
show_full(norm_constr_prob)
show_full(empirical_T)

latexify(round.(norm_constr_prob, digits=3))


## construct transition matrix directly
norm_constr_T, empirical_prob = construct_t(number_of_states, state_list[2500:end], holding_times[2500:end])
show_full(norm_constr_T)
show_full(round.(empirical_prob, digits=3))

latexify(round.(norm_constr_T, digits=3))
latexify(round.(empirical_prob, digits=3))
##

ss = steady_state(norm_constr_T, number_of_states=number_of_states)

round.(real.(ss), digits=3)

latexify(round.(real.(ss), digits=3))

### spectral bisection
#1 construct adjacency matrix
function get_distance_matrix(points)
    size = length(points)
    distance_matrix = Array{Any, 2}(nothing, (size, size))
    for i in range(1, size)
        for j in range(i, size)
            distance_matrix[i,j] = LinearAlgebra.norm(points[i]-points[j])
            distance_matrix[j,i] = LinearAlgebra.norm(points[i]-points[j])
        end
    end
    return distance_matrix
end

# get skipped points from run
function get_sample_points(X, Y, Z; step=100)
    points = []
    # size = length(X)
    # print(size)
    for i in range(start=10, stop=30000, step=step) #cut off beginning errata
        #print(i)
        #push!(points,[X[Int(i)],Y[Int(i)],Z[Int(i)]])
        push!(points,[X[i],Y[i],Z[i]])
    end
    return points
end

function get_graph_lap(distance_matrix; thresh=0.5)
    thresh = quantile(distance_matrix[:], thresh) #use the median as a threshold
    adj_mat = I - (distance_matrix .< thresh) #line borrowed from andre
    graph_lap = adj_mat - Diagonal(sum(adj_mat, dims = 1)[:])
    return graph_lap
end

sample_points = get_sample_points(X, Y, Z; step=100)
distance_matrix = get_distance_matrix(sample_points)
graph_lap = get_graph_lap(distance_matrix; thresh=0.7)

Λ, V = eigen(graph_lap)

function get_bisection_indices(V, index; thresh=0)
    bisection_indices = V[:, index] .> thresh
    return bisection_indices
end

function get_partition(V, points; thresh=0)
    #for now, just the one where we use 2,3,4
    list_2 = get_bisection_indices(V, 2, thresh=thresh)
    list_3 = get_bisection_indices(V, 3, thresh=thresh)
    list_4 = get_bisection_indices(V, 4, thresh=thresh)

    cluster_dict = Dict()
    map = Dict([0,0,0] => 1, [1,0,0] => 2, [0,1,0] => 3, [1,1,0] => 4,
                [0,0,1] => 5, [1,0,1] => 6, [0,1,1] => 7, [1,1,1] => 8)

    for i in range(1,length(points))
        cluster_number = map[[list_2[i], list_3[i], list_4[i]]]
        try
            cluster_dict[cluster_number] = [cluster_dict[cluster_number] points[i]]
        catch 
            cluster_dict[cluster_number] = points[i]
        end
    end
    #check to make sure all states exist

    # show_full(cluster_dict[1])
    # show_full(mean(cluster_dict[1], dims=2))
    out = Dict()
    for i in range(1,8) #once again, not general!
        try
            elements = cluster_dict[i]
            out[i] = vec(mean(elements,dims=2)) #centroid location
        catch KeyError
        end
        
    end
    return out
end

partition_dict = get_partition(V, sample_points)
centroids = [value for (key, value) in partition_dict]

Plots.scatter(Tuple.(sample_points), markersize=2)
Plots.scatter!(Tuple.(centroids))

# 
k_means = Clustering.kmeans(points, 9)
k_means_12 = Clustering.kmeans(points, 12)

Plots.scatter(Tuple.(sample_points), markersize=2, markerstrokewidth=0.5)
Plots.scatter!(Tuple.([x for x in eachcol(k_means_12.centers)]))

##judge different partitions

function create_thresh(partition_dict)
    # return a thresholding function that uses nearest-centroid from a dictionary of centroids
    function thresh_func(x,y,z)
        partition_dict = partition_dict
        min_dist = Inf64
        min_i = 0
        for i in range(1,8)
            try
                dist = LinearAlgebra.norm(partition_dict[i]-[x,y,z])
                if dist < min_dist
                    min_dist = dist
                    min_i = i
                end
            catch
            end
        end
        return min_i
    end
    return thresh_func
end

function test_partition(thresh_func)
    # run a simulation with that threshold functions
    X, Y, Z, state_list, holding_times = run_sim(;runs=1000000, timing=true, 
                            thresh_func=thresh_func, delta_rho=0, rho_start=28)
    # get T matrix
    norm_constr_T, empirical_prob = construct_t(number_of_states, state_list, holding_times)
    # return steady state solution
    return steady_state(norm_constr_T)
end


spectral_func = create_thresh(partition_dict)






## identify quantiles in Z
z_sort = copy(Z)
three_state_extremes = quantile!(z_sort,[0.05, 0.95]) #output: [10.6, 37.9]
five_state_extremes = quantile!(z_sort,[0.05, 0.25, 0.75, 0.95]) #output: [10.6, 16.9, 30.6, 37.9]
median = quantile!(z_sort,0.5)

###################

## basic plotting
using Plots
Plots.plot(X,Y)
Plots.plot(Z[100:end])
Plots.histogram(Z;bins=100)
##
