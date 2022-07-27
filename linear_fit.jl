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

elem_dict = get_elem_dict(T_dict)

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

lin_fit_matrix = get_lin_fit_matrix(elem_dict; number_of_states=12)

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