using Test
include("generalfunctions.jl")
# Definimos la ecuación diferencial
f(t, x, c) = x

# Parámetros para sol_ODE_param
c = ones(100) # puedes ajustar esto según tus necesidades
x₀ = 1.0
t₀ = 0.0
h = 0.01
N = 100
method = "RK4"

# Llamamos a sol_ODE_param
x, t = sol_ODE_param(f, c, x₀, t₀, h, N, method = method)

# Comprobamos que la solución es correcta
# La solución analítica para esta ecuación diferencial es e^t
@test all(abs.(x .- exp.(t)) .< 1e-3)

