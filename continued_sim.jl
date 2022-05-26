using Statistics

function show_full(array)
    show(IOContext(stdout, :limit=>false), MIME"text/plain"(), array)
end

## basic simulation functions
function xdot(x,y,z; sigma=10)
    return sigma*(y-x)
end
function ydot(x,y,z; rho=28)
    return x*(rho-z)-y
end
function zdot(x,y,z; beta=8/3)
    return x*y - beta*z
end

function make_step(x,y,z; delta=0.01)
    xp = x + xdot(x,y,z)*delta
    yp = y + ydot(x,y,z)*delta
    zp = z + zdot(x,y,z)*delta

    xn = x + 0.5*(xdot(x,y,z)+xdot(xp,yp,zp))*delta
    yn = y + 0.5*(ydot(x,y,z)+ydot(xp,yp,zp))*delta
    zn = z + 0.5*(zdot(x,y,z)+zdot(xp,yp,zp))*delta

    return xn, yn, zn
end

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
    


## simulation run
function run_sim(;runs=3000, timing=false, thresh_func=z_tercile_thresh)
    x, y, z = 0, 1, 0
    X, Y, Z = Float64[0,], Float64[1,], Float64[0,]
    cnt = 0
    state_list = [] #list of unique states
    holding_times = [] #list of holding times of each unique state

    state = thresh_func(x,y,z)

    for i in 1:runs
        x,y,z = make_step(x,y,z)
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

    end
    return X, Y, Z, state_list, holding_times
end

X, Y, Z, state_list, holding_times = run_sim(;runs=1000000, timing=true, thresh_func=twelve_state_quant_thresh)
number_of_states = 12

#basic plotting
using Plots
Plots.plot(X,Y)
Plots.plot(Z[100:end])
Plots.histogram(Z;bins=100)


#construct probability transition matrix directly
function construct_p(number_of_states, state_list, holding_times)
    constr_prob = zeros(number_of_states, number_of_states)

    for i in 1:length(state_list)-1
        local current = state_list[i]
        local next = state_list[i+1]
        constr_prob[current, current] += holding_times[i]-1
        constr_prob[next, current] += 1 #reconstruct full markov chain
    end
    norm = sum(constr_prob, dims=1) 

    norm_constr_prob = constr_prob ./ norm
    empirical_T = log(norm_constr_prob)/0.01
    return (norm_constr_prob, empirical_T)
end
norm_constr_prob, empirical_T = construct_p(number_of_states, state_list, holding_times)
show_full(norm_constr_prob)
show_full(empirical_T)

## other method! 
function construct_t(number_of_states, state_list, holding_times)
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

    empirical_prob = exp(norm_constr_T * 0.01)
    return (norm_constr_T, empirical_prob)
end

norm_constr_T, empirical_prob = construct_t(number_of_states, state_list, holding_times)
show_full(norm_constr_T)
show_full(empirical_prob)

`                                                           `

## identify quantiles in Z
z_sort = copy(Z)
three_state_extremes = quantile!(z_sort,[0.05, 0.95]) #output: [10.6, 37.9]
five_state_extremes = quantile!(z_sort,[0.05, 0.25, 0.75, 0.95]) #output: [10.6, 16.9, 30.6, 37.9]


###################

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
