
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
    if z<10.6 #bottom 
        return 1
    end
    if z<37.9 #middle
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
function twelve_state_just_high(x,y,z; mid=22.5, high=37.9)
    #using 0.5, 0.95 quantiles for z
    state = [0,0,0]
    if x > 0
        state[1] = 1
    end
    if y > 0 
        state[2] = 1
    end
    if z > high
        state[3] = 2
    elseif z > mid
        state[3] = 1
    end
    # inefficient to be defining the map inside the function...
    map = Dict([0,0,0] => 1, [1,0,0] => 2, [0,1,0] => 3, [1,1,0] => 4,
                [0,0,1] => 5, [1,0,1] => 6, [0,1,1] => 7, [1,1,1] => 8,
                [0,0,2] => 9, [1,0,2] => 10, [0,1,2] => 11, [1,1,2] => 12)
    return map[state]
end

function twelve_state_26(x,y,z)
    #using 0.5, 0.95 quantiles for z
    state = [0,0,0]
    if x > 0
        state[1] = 1
    end
    if y > 0 
        state[2] = 1
    end
    if z > 35.4
        state[3] = 2
    elseif z > 20.4
        state[3] = 1
    end
    # inefficient to be defining the map inside the function...
    map = Dict([0,0,0] => 1, [1,0,0] => 2, [0,1,0] => 3, [1,1,0] => 4,
                [0,0,1] => 5, [1,0,1] => 6, [0,1,1] => 7, [1,1,1] => 8,
                [0,0,2] => 9, [1,0,2] => 10, [0,1,2] => 11, [1,1,2] => 12)
    return map[state]
end

function spectral_thresh(x,y,z)
    centroids = Dict(      5 => [15.7328, 16.7556, 35.9377],
    4 => [3.43148, 2.73087, 24.7724],
    6 => [-5.51011, -5.539, 22.685],
    7 => [14.7184, 8.44782, 40.7711],
    2 => [2.94363, 5.61358, 17.6816],
    8 => [1.4769, -2.17442, 30.9628],
    3 => [10.1125, 9.07715, 29.8761],
    1 => [8.07909, 13.3689, 17.2107])
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
        