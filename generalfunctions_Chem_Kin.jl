"""
    sol_diff_eq_ord(f, x₀, t₀, h, N; method="RK4")

Resuelve una ecuación diferencial ordinaria utilizando varios métodos.

# Argumentos
- `f::Function`: Función de las variables t y x a evaluar.
- `x₀::Real`: Valor inicial de x.
- `t₀::Real`: Valor inicial de t.
- `h::Real`: Elemento diferencial dt.
- `N::Integer`: Número de pasos de t.
- `method::String="RK4"`: Método utilizado para evaluar la EDO.

# Devoluciones
- `x::Vector{Real}`: Vector de soluciones x.
- `t::Vector{Real}`: Vector de tiempos t.
"""
function sol_diff_eq_ord(f::Function, x₀::Real, t₀::Real, h::Real, N::Integer; method::String="RK4")
    t = range(t₀, step=h, length=N+1)
    x = zeros(Real, length(t))
    x[1] = x₀

    for n in 2:length(t)
        if method == "Euler"
            x[n] = x[n-1] + h*f(t[n-1], x[n-1])
        elseif method == "MidPoint"
            x[n] = x[n-1] + h*f(t[n-1] + 0.5*h, x[n-1] + 0.5*h*f(t[n-1], x[n-1]))
        elseif method == "Ralston"
            k₁ = f(t[n-1], x[n-1])
            k₂ = f(t[n-1] + (2/3)*h, x[n-1] + (2/3)*h*k₁)
            x[n] = x[n-1] + h*(0.25*k₁ + 0.75*k₂)
        elseif method == "RK4"
            k₁ = f(t[n-1], x[n-1])
            k₂ = f(t[n-1] + 0.5*h, x[n-1] + 0.5*h*k₁)
            k₃ = f(t[n-1] + 0.5*h, x[n-1] + 0.5*h*k₂)
            k₄ = f(t[n-1] + h, x[n-1] + h*k₃)
            x[n] = x[n-1] + h*(k₁ + 2*k₂ + 2*k₃ + k₄)/6
        else
            error("Método no reconocido: $method")
        end
    end

    return x, t
end
