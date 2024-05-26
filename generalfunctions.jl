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

"""
    sol_diff_eq_ord(f₁, f₂, x₁₀, x₂₀, t₀, h, N; method="RK4")

Resuelve dos ecuaciones diferenciales ordinarias simultáneamente.

# Argumentos
- `f₁`: Función que define la primera ecuación diferencial.
- `f₂`: Función que define la segunda ecuación diferencial.
- `x₁₀`: Valor inicial para la primera ecuación.
- `x₂₀`: Valor inicial para la segunda ecuación.
- `t₀`: Tiempo inicial.
- `h`: Paso de tiempo.
- `N`: Número de pasos.
- `method`: Método para resolver las ecuaciones (por defecto "RK4").

# Devoluciones
- `x₁`: Solución de la primera ecuación diferencial.
- `x₂`: Solución de la segunda ecuación diferencial.
- `t`: Puntos de tiempo correspondientes.
"""
function sol_two_diff_eq_ord(f₁::Function, f₂::Function, x₁₀::Real, x₂₀::Real, t₀::Real, h::Real, N::Integer; method::String="RK4")
    t = range(t₀, step=h, length=N+1)
    x₁ = zeros(Real, length(t))
    x₂ = zeros(Real, length(t))
    x₁[1] = x₁₀
    x₂[1] = x₂₀

    for n in 2:length(t)
        if method == "Euler"
            x₁[n] = x₁[n-1] + h*f₁(t[n-1], x₁[n-1])
            x₂[n] = x₂[n-1] + h*f₂(t[n-1], x₂[n-1])
        elseif method == "MidPoint"
            x₁[n] = x₁[n-1] + h*f₁(t[n-1] + 0.5*h, x₁[n-1] + 0.5*h*f₁(t[n-1], x₁[n-1]))
            x₂[n] = x₂[n-1] + h*f₂(t[n-1] + 0.5*h, x₂[n-1] + 0.5*h*f₂(t[n-1], x₂[n-1]))
        elseif method == "Ralston"
            k₁₁ = f₁(t[n-1], x₁[n-1])
            k₁₂ = f₁(t[n-1] + (2/3)*h, x₁[n-1] + (2/3)*h*k₁₁)
            x₁[n] = x₁[n-1] + h*(0.25*k₁₁ + 0.75*k₁₂)
            k₂₁ = f₂(t[n-1], x₂[n-1])
            k₂₂ = f₂(t[n-1] + (2/3)*h, x₂[n-1] + (2/3)*h*k₂₁)
            x₂[n] = x₂[n-1] + h*(0.25*k₂₁ + 0.75*k₂₂)
        elseif method == "RK4"
            k₁₁ = f₁(t[n-1], x₁[n-1])
            k₁₂ = f₁(t[n-1] + 0.5*h, x₁[n-1] + 0.5*h*k₁₁)
            k₁₃ = f₁(t[n-1] + 0.5*h, x₁[n-1] + 0.5*h*k₁₂)
            k₁₄ = f₁(t[n-1] + h, x₁[n-1] + h*k₁₃)
            x₁[n] = x₁[n-1] + h*(k₁₁ + 2*k₁₂ + 2*k₁₃ + k₁₄)/6
            k₂₁ = f₂(t[n-1], x₂[n-1])
            k₂₂ = f₂(t[n-1] + 0.5*h, x₂[n-1] + 0.5*h*k₂₁)
            k₂₃ = f₂(t[n-1] + 0.5*h, x₂[n-1] + 0.5*h*k₂₂)
            k₂₄ = f₂(t[n-1] + h, x₂[n-1] + h*k₂₃)
            x₂[n] = x₂[n-1] + h*(k₂₁ + 2*k₂₂ + 2*k₂₃ + k₂₄)/6
        else
            error("Método no reconocido: $method")
        end
    end

    return x₁, x₂, t
end



