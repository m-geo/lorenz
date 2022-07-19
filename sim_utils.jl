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