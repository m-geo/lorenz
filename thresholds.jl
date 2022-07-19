
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

function spectral_thresh(x,y,z)
    centroids = Dict(  5 => [9.98259, 13.3079, 23.0984],
                        4 => [-11.6608, -12.1595, 29.5316],
                        6 => [3.54258, 5.51565, 16.788],
                        7 => [14.288, 11.9913, 36.8665],
                        2 => [-3.35442, -3.47187, 22.2316],
                        8 => [-17.0136, -20.4666, 35.5536],
                        3 => [12.6386, 4.59547, 39.5217],
                        1 => [7.58452, 0.825746, 33.1679])
    min_dist = Inf64
    min_i = 0
    for i in range(1,8)
        dist = LinearAlgebra.norm(centroids[i]-[x,y,z])
        if dist < min_dist
            min_dist = dist
            min_i = i
        end
    end
    return min_i
end
        