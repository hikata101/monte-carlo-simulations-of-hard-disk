def periodic(condition, value, debug):
    original_condition = condition  # Preserve the original condition value
    if value > original_condition:
        if debug:
            print("subtracted", value, "by", original_condition)
        value = value - original_condition * 2
        if debug:
            print(value)
    if value < -original_condition:
        if debug:
            print("added", value, "by", original_condition)
        value = value + (original_condition * 2)
        if debug:
            print(value)
    return value

def x_evolve_pbc(gamma, dt, g, former_x, length, debug):
    next_x = former_x+math.sqrt(2*gamma*dt)*g
    return periodic(length, next_x, debug)

def x_list_pbc(gamma, dt, n, length, debug):
    condition = length / 2
    tmp_x = np.random.normal(0, 33.4, 1)
    tmp_x = periodic(tmp_x, condition, debug)
    result = [tmp_x]
    for i in np.random.normal(0, 1, n):
        tmp_x = x_evolve_pbc(1,1,i,tmp_x, condition, debug)
        result.append(tmp_x)
    return result