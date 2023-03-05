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

function show_full(array)
    show(IOContext(stdout, :limit=>false), MIME"text/plain"(), array)
end

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

##

function new_run_sim(;runs=3000, timing=false, thresh_func=z_tercile_thresh, delta_rho=0, rho_start=28)
    x, y, z = 14.0, 20.0, 27.0
    X, Y, Z = Float64[0,], Float64[1,], Float64[0,]
    cnt = 0
    state_list = [] #list of unique states
    holding_times = [] #list of holding times of each unique state

    sim_list = []  #list of succesive 3-tuples of the system's location
    markov_chain = [] #simple list of successive states

    state = thresh_func(x,y,z)
    rho = rho_start
    rho_step = delta_rho / runs

    for i in 1:runs
        #print(rho)
        x,y,z = make_step(x,y,z; rho=rho)
        push!(X, x)
        push!(Y, y)
        push!(Z, z)
        state = [x,y,z]
        push!(sim_list, state)
        push!(markov_chain, thresh_func(x,y,z))

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
    # return X, Y, Z, state_list, holding_times
    return sim_list, markov_chain
end


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