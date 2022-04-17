using Statistics


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
function state_function(x,y,z)
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

function threshold_basic(x,y,z)
    # Z=24
    out = 0
    if z>24
        out = 1
    end
    return out
end

## simulation run
function run_sim(;runs=3000, timing=false, thresh_func=state_function)
    x, y, z = 0, 1, 0
    X, Y, Z = Float64[0,], Float64[1,], Float64[0,]
    cnt = 0
    #cnt_list = Tuple[]
    state_list = []
    holding_times = []

    state = thresh_func(x,y,z)

    for i in 1:runs
        x,y,z = make_step(x,y,z)
        push!(X, x)
        push!(Y, y)
        push!(Z, z)

        if thresh_func(x,y,z) == state
            cnt += 1
        else
            #push!(cnt_list, (state, cnt, thresh_func(x,y,z))) #old state, time held, new state
            push!(state_list, state)
            push!(holding_times, cnt)
            state = thresh_func(x,y,z)
            cnt = 0
        end

    end
    return X, Y, Z, state_list, holding_times
end

X, Y, Z, state_list, holding_times = run_sim(;runs=300000, timing=true, thresh_func=state_function)
number_of_states = 3

#basic plotting
using Plots
Plots.plot(X,Y)
Plots.plot(Z[100:end])
Plots.histogram(Z;bins=100)

#times is a list of 3-elem tuples
#Plots.histogram([x[2]*0.01 for x in times[3:end] if x[1]==0],bins=20)

#construct probability transition matrix directly
constr_prob = zeros(number_of_states, number_of_states)

for i in 1:length(state_list)-1
    local current = state_list[i]
    local next = state_list[i+1]
    constr_prob[current, current] += holding_times[i]-1
    constr_prob[next, current] += 1 #reconstruct full markov chain
end

norm = sum(constr_prob, dims=1) 
norm_constr_prob = constr_prob ./ norm
# --> okay so this makes sense

## other method! 
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
    norm_constr_T[i, i] = -1.0
    norm_constr_T[:, i] *= holding_scale[i]
end
empirical_constr_T = norm_constr_T #think about what this is

empirical_prob = exp(empirical_constr_T * 0.01)

# so this result feels correct ^^ (the other one doesn't)



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
