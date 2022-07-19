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
