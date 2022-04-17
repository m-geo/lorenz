
function xdot(x,y,z; sigma=10)
    return sigma*(y-x)
end

function ydot(x,y,z; rho=48)
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



function run_sim(;runs=3000, count_z = false, timing=false, thresh = 0.0)
    x, y, z = 0, 1, 0
    X, Y, Z = Float64[0,], Float64[1,], Float64[0,]
    #dist = fill(0,2)
    cnt = 0
    zless = true
    cnt_list = Int[]

    for i in 1:runs
        x,y,z = make_step(x,y,z)
        push!(X, x)
        push!(Y, y)
        push!(Z, z)

        if (z < thresh && zless) || (z > thresh && !zless)
            cnt += 1
        else
            push!(cnt_list, cnt)
            zless = !zless
            cnt = 0
        end

        # if count_z
        #     if z > thresh
        #         dist[1] += 1
        #     else
        #         dist[2] += 1
        #     end
        # end
    end
    return X, Y, Z, cnt_list #,dist
end

X, Y, Z, times = run_sim(;runs=30000, timing=true, count_z = false, thresh = 24)

using Plots
Plots.plot(X,Y)
Plots.plot(Z)
Plots.histogram(Z;bins=100)

Plots.histogram((times*0.01)[3:end],bins=20)

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

# get subsequent maxima plot
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

#=
points = Observable(Point3f[])
colors = Observable(Int[])
set_theme!(theme_black())

fig, ax, l = lines(points, color = colors,
    colormap = :inferno, transparency = true,
    axis = (; type = Axis3, protrusions = (0, 0, 0, 0),
        viewmode = :fit, limits = (-30, 30, -30, 30, 0, 50)))


frames = 1:120
for frame in frames
    for i in 1:50
        push!(points[], step!(attractor))
        push!(colors[], frame)
    end
    ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120)
    notify.((points, colors))
    l.colorrange = (0, frame)
    sleep(0.1)
    #println(frame)
end
=#
