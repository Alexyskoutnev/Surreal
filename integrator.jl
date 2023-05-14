using Plots

function f_exact(t, x)
    return sin(t) + x^2 / 2
end

function f(t::Float64, x::Float64)
    """Example derivative function x' = sin(y) + x"""
    return cos(t) + x
end

function forward_euler(x, t, step_size)
    for i in 1:num-1
        x[i + 1] = x[i] + step_size * f(t[i], x[i])
    end
    return t, x
end

function main()
    t_start = 0
    num = 10000
    t_end = 100
    t = range(t_start,stop=t_end,length=num)
    step_size = (t_start - t_end) / num 
    x_0, t_0 = 1.0, t[1]
    x = fill(1.0, num)
    x_real = fill(1.0, num)
    #Initial Condition
    x[1] = x_0
    x_real[1] = x_0
    #Numerical Approximation
    t, x_approx = forward_euler(x, t, step_size)
    i = 1
    for (i, _t) in enumerate(t[2:end])
        val = f_exact(_t, x_real[i])
        x_real[i + 1] = val
    end
    plot(t, x, xlabel="time", ylabel="f(x)")
end