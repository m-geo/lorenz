using Polynomials

### explore fitting function to T matrix
function get_T_dict(;runs=1000000, thresh_func=twelve_state_just_high, number_of_states=12, 
                    rho_range=range(28,step=4,stop=48))
    # generates a disctionary keyed on rho values and valued with the associated transiton matrivces
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

T_dict = get_T_dict(runs=1000000, thresh_func=twelve_state_just_high, number_of_states=12, 
                    rho_range=range(28,step=2,stop=48))

Plots.heatmap(T_dict[48], yflip=true, clim=(-80,80))

latexify(round.(T_dict[32], digits=3))


function get_elem_dict(T_dict; number_of_states=12, rho_range=range(28,step=4,stop=48))
    #this is NOT the most elegant way to do it 
    #better treat it not as a dict but as a 3d array so yuo can index it directly
    # divides the T_dict into a dictionary of index number: list of entries at that index
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

elem_dict = get_elem_dict(static_generators; rho_range=range(26,step=1,stop=32))

function get_lin_fit_matrix(elem_dict; number_of_states=12, rho_range=range(28,step=4,stop=48))
    # using the elementwise dictionary, generate a matrix of linear fit polynomials
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

lin_fit_matrix = get_lin_fit_matrix(elem_dict; number_of_states=12,  rho_range=range(26,step=1,stop=32))

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
function test_fit(lin_fit_matrix; rho_range=range(start=29,step=2,stop=47), number_of_states=12)
    distance_array = []
    for rho in rho_range
        # println(rho)
        X, Y, Z, state_list, holding_times = run_sim(;runs=1000000, timing=true, 
                    thresh_func=twelve_state_just_high, delta_rho=0, rho_start=rho)
        # println("sim done!")
        norm_constr_T, empirical_prob = construct_t(number_of_states, state_list, holding_times)
        fitted_matrix = fit_matrix(lin_fit_matrix, rho; number_of_states)
        distance = LinearAlgebra.norm(norm_constr_T-fitted_matrix)/LinearAlgebra.norm(norm_constr_T)
        # println(distance)
        push!(distance_array, distance)
    end
    return distance_array
end

distance_array = test_fit(lin_fit_matrix)       # mean is 0.33


# ---------- new ---------
static_generators = JLD.load("./vars/static_generators.jld", "static_generators")
further_generators = JLD.load("./vars/further_generators.jld", "further_generators")

static_means = JLD.load("./vars/bayesian_static_means.jld",  "bayesian_static_means")
static_stds = JLD.load("./vars/bayesian_static_stds.jld",  "bayesian_static_stds")

function plot_diagonals(generator_dict, static_means, static_stds; extension_dict=nothing, sliding_windows=nothing, middle_values=nothing, bayesian_delta=nothing)
    ref_list = [x for x in 26:32]
    if extension_dict !== nothing
        max_rho = 36

        full_ref_list = [x for x in 26:36]

    else
        max_rho=32
    end


    fig = Figure(resolution=(1600,1200))
    for i in 1:3, j in 1:4
        box = j+4*(i-1) 
        ax = Axis(fig[i,j], title="$box")
        entry_list = [generator_dict[rho][box,box] for rho in ref_list]

        mean_list = [] # this is sooo not general
        std_list = []
        for rho in 26:1:max_rho
            push!(mean_list, static_means[rho][box, box])
            push!(std_list, static_stds[rho][box, box])
        end

        f1 = [Polynomials.fit(ref_list, entry_list, 1)[i] for i in 0:1]
        f2 = [Polynomials.fit(ref_list, entry_list, 2)[i] for i in 0:2]

        if extension_dict !== nothing
            for rho in 33:36
                push!(entry_list, extension_dict[rho][box,box])
            end
            scatter!(full_ref_list, entry_list)
            errorbars!(full_ref_list, mean_list, std_list)
            xs = full_ref_list[1]:0.1:full_ref_list[end]
            lines!(xs, [f1[1]+f1[2]*x for x in xs])
            lines!(xs, [f2[1]+f2[2]*x+f2[3]*x^2 for x in xs])
        else
            scatter!(ref_list, entry_list)
        end

        if sliding_windows !== nothing
            sliding_means = [mean(sliding_windows[i])[box, box] for i in 1:1:5]
            sliding_stds = [std(sliding_windows[i])[box, box] for i in 1:1:5]
            errorbars!(middle_values, sliding_means, sliding_stds, color="red")
            scatter!(middle_values, sliding_means, color="red")
        end

        if bayesian_delta !== nothing
            scatter!(29, mean(bayesian_delta)[box,box],color="blue")
            errorbars!([29], [mean(bayesian_delta)[box,box]], [std(bayesian_delta)[box,box]], color="blue")
        end

    end
    fig
end

#inputs here taken from smoothness_test.jl

sim_list_delta, markov_chain_delta = new_run_sim(;runs=1e6, timing=true, 
                            thresh_func=twelve_state_just_high, delta_rho=6, rho_start=26)

bayesian_delta = BayesianGenerator(markov_chain_delta; dt=dt)

plot_diagonals(static_generators, static_means, static_stds; #extension_dict=further_generators, 
    sliding_windows=sliding_bayesians, middle_values=middle_values, bayesian_delta=bayesian_delta)

function plot_transitions(generator_dict)
    ref_list = [x for x in 26:32]
    fig = Figure(resolution=(800,400))
    box_list = [(5,9), (8,12)]
    for n in eachindex(box_list)
        i, j = box_list[n]
        ax = Axis(fig[1,n], title="$i -> $j")
        entry_list = [generator_dict[rho][j,i] for rho in ref_list] #note j and i flipped
        scatter!(ref_list, entry_list)
    end
    fig
end

plot_transitions(static_generators)