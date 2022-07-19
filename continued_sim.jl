using Statistics
using DataFrames
using Clustering
using Polynomials
using LinearAlgebra
using Makie
include("thresholds.jl")
include("sim_utils.jl")

function show_full(array)
    show(IOContext(stdout, :limit=>false), MIME"text/plain"(), array)
end

#run simulation:
X, Y, Z, state_list, holding_times = run_sim(;runs=1000000, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=20, rho_start=28)
number_of_states = 8
##

## construct probability matrix directly
function construct_p(number_of_states, state_list, holding_times; simple=false)
    constr_prob = zeros(number_of_states, number_of_states)

    for i in 1:length(state_list)-1
        local current = state_list[i]
        local next = state_list[i+1]
        constr_prob[current, current] += holding_times[i]-1
        constr_prob[next, current] += 1 #reconstruct full markov chain
    end
    norm = sum(constr_prob, dims=1) 
    show_full(constr_prob)
    show_full(norm)
    norm_constr_prob = constr_prob ./ norm
    
    #set columns with negative values to zero ()
    for col in eachcol(norm_constr_prob)
        #print(col)
        for x in col
            if isnan(x)
                fill!(col, 0.0)
                break
            end
        end
    end

    if simple
        return (norm_constr_prob, NaN)
    else
        # # try
        # #     empirical_T = log(norm_constr_prob)/0.01
        # # catch e
        # empirical_T =  copy(norm_constr_prob)
        # #print("we're in the catch statement" + empirical_T)
        # for col in eachcol(norm_constr_prob)
        #     if !all(x->x==0.0,col)
        #         println("found!")
        #         print(col)
        #         col = log(col)/0.01
        #     end
        # end
        # # end
        empirical_T = log(norm_constr_prob)/0.01
        return (norm_constr_prob, empirical_T)
    end
end
norm_constr_prob, empirical_T = construct_p(number_of_states, state_list, holding_times, simple=false)
show_full(norm_constr_prob)
show_full(empirical_T)

latexify(round.(norm_constr_prob, digits=3))

## construct transition matrix directly
function construct_t(number_of_states, state_list, holding_times; simple=false)
    holding_time_dist = [[] for n in 1:number_of_states]
    constr_T = zeros(number_of_states, number_of_states)

    #count number of times unique states appear
    for i in 1:length(state_list)-1 #for i in unique states #why is the -1??
        local current = state_list[i]
        local next = state_list[i+1]
        constr_T[next, current] += 1 
        push!(holding_time_dist[current], holding_times[i]*0.01)
    end

    norm = sum(constr_T, dims=1)
    norm_constr_T = constr_T ./ norm

    holding_scale = 1 ./ mean.(holding_time_dist) 
    for i in 1:number_of_states
        norm_constr_T[i, i] = -1.0 #this is the bit I don't quite understand
        norm_constr_T[:, i] *= holding_scale[i]
    end

    if simple
        return norm_constr_T
    else
        empirical_prob = exp(norm_constr_T * 0.01)
        return (norm_constr_T, empirical_prob)
    end
end
norm_constr_T, empirical_prob = construct_t(number_of_states, state_list, holding_times)
show_full(norm_constr_T)
show_full(empirical_prob)

latexify(round.(norm_constr_T, digits=3))
##


### explore fitting function to T matrix
function get_T_dict(;runs=1000000, thresh_func=twelve_state_just_high, number_of_states=12, rho_range=range(28,step=4,stop=48))
    T_dict = Dict()
    for rho in rho_range
        X, Y, Z, state_list, holding_times = run_sim(;runs=runs, timing=true, 
                                            thresh_func=thresh_func, delta_rho=0, rho_start=rho)
        norm_constr_T, empirical_prob = construct_t(number_of_states, state_list, holding_times)
        #push!(T_list, (rho,norm_constr_T))
        T_dict[rho] = norm_constr_T
    end
    return T_dict
end

T_dict = get_T_dict()

Plots.heatmap(T_dict[48], yflip=true, clim=(-80,80))

latexify(round.(T_dict[32], digits=3))


function get_elem_dict(T_dict; number_of_states=12, rho_range=range(28,step=4,stop=48))
    #this is NOT the most elegant way to do it 
    #better treat it not as a dict but as a 3d array so yuo can index it directly
    elem_dict = Dict()
    for i in range(1,number_of_states^2)
        for rho in rho_range
            try
                push!(elem_dict[i], T_dict[rho][i])
            catch
                elem_dict[i] = [T_dict[rho][i]]
            end
        end
    end
    return elem_dict
end

elem_dict = get_elem_dict(T_dict)

function get_lin_fit_matrix(elem_dict; number_of_states=12, rho_range=range(28,step=4,stop=48))
    lin_fit_matrix = Array{Any,2}(nothing,(number_of_states,number_of_states))
    rho_list = [i for i in rho_range]
    for i in range(1,number_of_states^2)
        lineq = Polynomials.polyfitA(rho_list, elem_dict[i],1)
        lin_fit_matrix[i] = lineq
    end
    for i in range(1, number_of_states)
        lin_fit_matrix[i,i] = 0.0
    end
    return lin_fit_matrix
end

lin_fit_matrix = get_lin_fit_matrix(elem_dict; number_of_states)

function fit_matrix(lin_fit_matrix, rho; number_of_states=12)
    t_matrix = zeros(number_of_states,number_of_states)
    #cnt = 0
    for i in range(1, number_of_states^2)
        #the entry is a function, so get func (rho)
        try
            t_matrix[i] = lin_fit_matrix[i](rho)
        catch e
            t_matrix[i] == 0.0
            # cnt += 1
        end
    end
    #print(cnt)
    for j in range(1, number_of_states)
        #set diagonal entries to the neg sum of column
        t_matrix[j,j] = -sum(t_matrix[:,j])
    end
    return t_matrix
end

# to test the lin fits

function test_fit(lin_fit_matrix; rho_range=(start=29,step=2,stop=47), number_of_states=12)
    distance_array = []
    for rho in rho_range
        print(rho)
        X, Y, Z, state_list, holding_times = run_sim(;runs=1000000, timing=true, 
                    thresh_func=twelve_state_just_high, delta_rho=0, rho_start=rho)
        norm_constr_T, empirical_prob = construct_t(number_of_states, state_list, holding_times)
        fitted_matrix = fit_matrix(lin_fit_matrix, rho; number_of_states)
        distance = LinearAlgebra.norm(norm_constr_T-fitted_matrix)/LinearAlgebra.norm(norm_constr_T)
        push!(distance_array, distance)
    end
    return mean(distance_array)
end

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
function get_sample_points(X, Y, Z)
    points = []
    # size = length(X)
    # print(size)
    for i in range(start=10, stop=30000, step=100) #cut off beginning errata
        #print(i)
        #push!(points,[X[Int(i)],Y[Int(i)],Z[Int(i)]])
        push!(points,[X[i],Y[i],Z[i]])
    end
    return points
end

function get_graph_lap(distance_matrix)
    thresh = quantile(distance_matrix[:], 0.5) #use the median as a threshold
    adj_mat = I - (distance_matrix .< thresh) #line borrowed from andre
    graph_lap = adj_mat - Diagonal(sum(adj_mat, dims = 1)[:])
    return graph_lap
end

sample_points = get_sample_points(X, Y, Z)
distance_matrix = get_distance_matrix(sample_points)
graph_lap = get_graph_lap(distance_matrix)

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
    # show_full(cluster_dict[1])
    # show_full(mean(cluster_dict[1], dims=2))
    out = Dict()
    for i in range(1,8) #once again, not general!
        elements = cluster_dict[i]
        out[i] = vec(mean(elements,dims=2)) #centroid location
    end
    return out
end

partition_dict = get_partition(V, sample_points)

#next step: plot the locations, relative to the butterfly

Plots.plot(Tuple.(sample_points))




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
