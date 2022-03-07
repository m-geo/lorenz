
function xdot(x,y,z; sigma=10)
    return sigma*(y-x)
end

function ydot(x,y,z; rho=28)
    return x*(rho-z)-y
end

function zdot(x,y,z; beta=8/3)
    return x*y - beta*z
end


function step(x,y,z; delta=0.01)
    xp = x + xdot(x,y,z)*delta
    yp = y + ydot(x,y,z)*delta
    zp = z + zdot(x,y,z)*delta

    xn = x + 0.5*(xdot(x,y,z)+xdot(xp,yp,zp))*delta
    yn = y + 0.5*(ydot(x,y,z)+ydot(xp,yp,zp))*delta
    zn = z + 0.5*(zdot(x,y,z)+zdot(xp,yp,zp))*delta

    return xn, yn, zn
end

# testing it out
x, y, z = 0, 1, 0
X, Y, Z = Float64[0,], Float64[1,], Float64[0,]

for i in 1:3000
    x,y,z = step(x,y,z)
    push!(X, x)
    push!(Y, y)
    push!(Z, z)
end

using Plots
Plots.plot(X,Y)
Plots.plot(Z)


function get_maxima(ls)
    maxs = Float64[]
    if ls[1]>ls[2]
        push!(maxs,ls[1])
    end
    for i=2:length(ls)-1
        if ls[i-1]<ls[i]>ls[i+1]
            push!(maxs,ls[i])
        end
    end
    if ls[end]>ls[end-1]
        push!(maxs,ls[end])
    end
    return maxs
  end

z_max = get_maxima(Z)
z_plot = Tuple{Float64, Float64}[]
for i=1:length(z_max)-1
    push!(z_plot,(z_max[i], z_max[i+1]))
end

i1 = [x[1] for x in z_plot]
i2 = [x[2] for x in z_plot]
Plots.plot(i1,i2)


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
