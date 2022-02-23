
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

    x = x + 0.5*(xdot(x,y,z)+xdot(xp,yp,zp))*delta
    y = y + 0.5*(ydot(x,y,z)+ydot(xp,yp,zp))*delta
    z = z + 0.5*(zdot(x,y,z)+zdot(xp,yp,zp))*delta

    return x,y,z
end

# testing it out
x, y, z = 0, 1, 0
X, Y, Z = Float64[0,], Float64[1,], Float64[0,]

for i in 1:10000
    x,y,z = step(x,y,z)
    push!(X, x)
    push!(Y, y)
    push!(Z, z)
end


