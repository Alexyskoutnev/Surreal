using Plots, Infiltrator


# function f_exact(t, x)
#     return sin(t) + x^2 / 2
# end

# function f(t::Float64, x::Float64)
#     """Example derivative function x' = sin(y) + x"""
#     return cos(t) + x
# end

f(x) = x^2
f_real(x) = x^3 / 3.0 + 1

function forward_euler(x::Float64; h = 0.001)
    x_k = x + h * f(x)
end

function forward_euler(x::Vector{Float64}; h = 0.001)
    x_k = x + h * f(x)
end

function simulate!(xtraj, N; h = 0.001)
    for k = 1:(N-1)
        xtraj[:, k+1] .= forward_euler(xtraj[k], h=h)
    end 
end

function main()
    T_f = 1
    h = 0.001
    N = Int(floor(T_f./h + 1))
    t = h.*Array(0:(N-1))

    #Numerical Approx
    x0 = 1
    xtraj = zeros(1, N)
    xtraj[1] = x0

    #Simulate/Solve the ODE
    simulate!(xtraj, N, h=h)
    xtraj_a = [f_real(x) for x in t]
    plot(t, vec(xtraj))
    plot!(t, xtraj_a)
end