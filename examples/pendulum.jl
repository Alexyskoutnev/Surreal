import Pkg; Pkg.add("Pkg")

using LinearAlgebra
using PyPlot
using Infiltrator

l = 1.0
m = 1.0
g = 9.81

function f(x)
    θ = x[1]
    θ̇ = x[2]
    ẋ = [θ̇; (-g / l)*sin(θ)]
end

function euler_step(xk)
    xn = xk + h * f(xk)
end

function simulate!(xtraj, N)
    for k = 1:(N-1)
        xtraj[:, k+1] .= euler_step(xtraj[:,k])
    end
end

Tf = 10.0
h = 0.05
N = Int(floor(Tf./h + 1))
thist = h.*Array(0:(N-1))

x0 = [30*(pi/180); 0.0]
xtraj = zeros(2, N)
xtraj[:,1] = x0


function main()
    simulate!(xtraj, N)
    plot(thist,xtraj[1,:])
    xlabel("Time (sec)")
    ylabel("θ (rad)")
end