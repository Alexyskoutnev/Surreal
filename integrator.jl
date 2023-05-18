using Plots, Infiltrator, ForwardDiff, LinearAlgebra

f(x) = x^2
f_real(x) = x^3 / 3.0 + 1

#Global constants
l = 1.0
m = 1.0
g = 9.81
c = 0.1

function f_damped_pendulum(x)
    θ = x[1]
    θ̇ = x[2]
    ẋ = [θ̇; (-g/l)*sin(θ)-c*θ̇]
    ẋ
end

function forward_euler(x::Float64; h = 0.001)
    x_k = x + h * f(x)
end

function forward_euler(x::Vector{Float64}; h = 0.001)
    x_k = x + h * f_damped_pendulum(x)
end

function taylor_step(x; h=0.001)
    Ak = ForwardDiff.jacobian(f_damped_pendulum, x)
    xk = x + h*f_damped_pendulum(x) + 0.5*h*h*Ak*f_damped_pendulum(x)
    xk
end

function explicit_midpoint_step(x; h=0.001)
    xm = x + 0.5*h*f_damped_pendulum(x)
    xn = x + h*f_damped_pendulum(xm)
    xn
end

function simulate_euler!(xtraj, N; h = 0.001)
    for k = 1:(N-1)
        xtraj[:, k+1] .= forward_euler(xtraj[:,k], h=h)
    end 
end

function simulate_taylor!(xtraj, N; h = 0.001)
    for k = 1:(N-1)
        xtraj[:, k+1] .= taylor_step(xtraj[:,k], h=h)
    end 
end

function simulate_explicit_midpoint_step!(xtraj, N; h = 0.001)
    for k = 1:(N-1)
        xtraj[:, k+1] .= explicit_midpoint_step(xtraj[:,k], h=h)
    end 
end

function test_ODE_Solvers()
    T_f = 10
    h = 0.01
    N = Int(floor(T_f./h + 1))
    t = h.*Array(0:(N-1))

    #Numerical Approx
    x0 =[45*(pi/180) ; 0.0]
    xtraj_euler = zeros(2, N)
    xtraj_taylor = zeros(2, N)
    xtraj_midpoint = zeros(2, N)
    xtraj_euler[:,1] = x0
    xtraj_taylor[:,1] = x0
    xtraj_midpoint[:,1] = x0

    #Solve ODE
    simulate_euler!(xtraj_euler, N, h=h)
    simulate_taylor!(xtraj_taylor, N, h=h)
    simulate_explicit_midpoint_step!(xtraj_midpoint, N, h=h)

    #Plotting the Damped Pendulum Trajectory
    plot(t, xtraj_euler[1,:], label="euler")
    plot!(t, xtraj_taylor[1,:], label="taylor")
    plot!(t, xtraj_midpoint[1,:], label="midpoint")
    xlabel!("time [s]")
    ylabel!("angle [θ]")

    #Angular Velocity
    # plot(t, xtraj_euler[2,:], label="euler")
    # plot!(t, xtraj_taylor[2,:], label="taylor")
    # plot!(t, xtraj_midpoint[2,:], label="midpoint")
    # xlabel!("time [s]")
    # ylabel!("angular velocity")

end
