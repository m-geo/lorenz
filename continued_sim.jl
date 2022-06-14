using Statistics
using DataFrames
using Clustering
using Polynomials

function show_full(array)
    show(IOContext(stdout, :limit=>false), MIME"text/plain"(), array)
end

## basic simulation functions
function xdot(x,y,z; sigma=10)
    return sigma*(y-x)
end
function ydot(x,y,z; rho=28)
    #print(rho)
    return x*(rho-z)-y  ##notice how rho-z is the avg temperature gradient of sorts
end
function zdot(x,y,z; beta=8/3)
    return x*y - beta*z
end
function make_step(x,y,z; rho = 28, delta=0.01)
    #print(rho)
    xp = x + xdot(x,y,z)*delta
    yp = y + ydot(x,y,z; rho)*delta
    zp = z + zdot(x,y,z)*delta

    xn = x + 0.5*(xdot(x,y,z)+xdot(xp,yp,zp))*delta
    yn = y + 0.5*(ydot(x,y,z; rho)+ydot(xp,yp,zp; rho))*delta
    zn = z + 0.5*(zdot(x,y,z)+zdot(xp,yp,zp))*delta

    return xn, yn, zn
end
##

## threshold functions
function threshold_basic(x,y,z)
    # Z=24
    out = 0
    if z>24
        out = 1
    end
    return out
end
function z_tercile_thresh(x,y,z)
    # states should be mutually exclusive
    # Z=18.3 and Z=32
    if z<18.3 #bottom 
        return 1
    end
    if z<32 #middle
        return 2
    end 
    return 3 #top
end
function eight_state_thresh(x,y,z)
    state = [0,0,0]
    if x > 0
        state[1] = 1
    end
    if y > 0 
        state[2] = 1
    end
    if z > 24
        state[3] = 1
    end
    map = Dict([0,0,0] => 1, [1,0,0] => 2, [0,1,0] => 3, [1,1,0] => 4,
                [0,0,1] => 5, [1,0,1] => 6, [0,1,1] => 7, [1,1,1] => 8)
    return  map[state]
end
function twelve_state_quant_thresh(x,y,z)
    #using 0.05, 0.95 quantiles for z
    state = [0,0,0]
    if x > 0
        state[1] = 1
    end
    if y > 0 
        state[2] = 1
    end
    if z > 37.9
        state[3] = 2
    elseif z > 10.6
        state[3] = 1
    end
    # inefficient to be defining the map inside the function...
    map = Dict([0,0,0] => 1, [1,0,0] => 2, [0,1,0] => 3, [1,1,0] => 4,
    [0,0,1] => 5, [1,0,1] => 6, [0,1,1] => 7, [1,1,1] => 8,
    [0,0,2] => 9, [1,0,2] => 10, [0,1,2] => 11, [1,1,2] => 12)
    return map[state]
end
function twelve_state_just_high(x,y,z)
    #using 0.5, 0.95 quantiles for z
    state = [0,0,0]
    if x > 0
        state[1] = 1
    end
    if y > 0 
        state[2] = 1
    end
    if z > 37.9
        state[3] = 2
    elseif z > 22.5
        state[3] = 1
    end
    # inefficient to be defining the map inside the function...
    map = Dict([0,0,0] => 1, [1,0,0] => 2, [0,1,0] => 3, [1,1,0] => 4,
                [0,0,1] => 5, [1,0,1] => 6, [0,1,1] => 7, [1,1,1] => 8,
                [0,0,2] => 9, [1,0,2] => 10, [0,1,2] => 11, [1,1,2] => 12)
    return map[state]
end
##  



## simulation run
function run_sim(;runs=3000, timing=false, thresh_func=z_tercile_thresh, delta_rho=0, rho_start=28)
    x, y, z = 0, 1, 0
    X, Y, Z = Float64[0,], Float64[1,], Float64[0,]
    cnt = 0
    state_list = [] #list of unique states
    holding_times = [] #list of holding times of each unique state

    state = thresh_func(x,y,z)
    rho = rho_start
    rho_step = delta_rho / runs

    for i in 1:runs
        #print(rho)
        x,y,z = make_step(x,y,z; rho=rho)
        push!(X, x)
        push!(Y, y)
        push!(Z, z)

        if thresh_func(x,y,z) == state
            cnt += 1
        else
            push!(state_list, state)
            push!(holding_times, cnt)
            state = thresh_func(x,y,z)
            cnt = 0
        end
        rho += rho_step
    end
    return X, Y, Z, state_list, holding_times
end

X, Y, Z, state_list, holding_times = run_sim(;runs=1000000, timing=true, thresh_func=twelve_state_just_high, delta_rho=0, rho_start=28)
number_of_states = 12
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


# list_of_first_elements = []
# for i in T_list
#     push!(list_of_first_elements, (i[1], i[2][1]))
# end
#Plots.scatter([x[1] for x in list_of_first_elements], [x[2] for x in list_of_first_elements])
function get_elem_dict(T_dict; number_of_states=12, rho_range=range(28,step=4,stop=48))
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
    return lin_fit_matrix
end

lin_fit_matrix = get_lin_fit_matrix(elem_dict; number_of_states)




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

function get_maxima(ls; add_min=false)
    maxs = Float64[]
    if add_min
        mins = Float64[]
        if ls[1]<ls[2]
            push!(mins,ls[1])
        end
    end
    if ls[1]>ls[2]
        push!(maxs,ls[1])
    end
    for i=2:length(ls)-1
        if ls[i-1]<ls[i]>ls[i+1]
            push!(maxs,ls[i])
        end
        if add_min
            if ls[i-1]>ls[i]<ls[i+1]
                push!(mins,ls[i])
            end
        end
    end
    if ls[end]>ls[end-1]
        push!(maxs,ls[end])
    end
    if add_min
        if ls[end]<ls[end-1]
            push!(mins,ls[end])
        end
        return maxs, mins
    else
        return maxs
    end
end

# get subsequent maxima  "tent" plot
z_max = get_maxima(Z)
z_plot = Tuple{Float64, Float64}[]
for i=1:length(z_max)-1
    push!(z_plot,(z_max[i], z_max[i+1]))
end
i1 = [x[1] for x in z_plot]
i2 = [x[2] for x in z_plot]
Plots.scatter(i1,i2)

# get sequential max/min plot
z_max, z_min = get_maxima(Z, add_min=true)
Plots.scatter(z_max, z_min)
