
include("generalfunctions.jl")
include("extrafunctions.jl")



"""
Dimerization(t₀::Real, χ₀::Vector{Float64}, k₁::Real, k₂::Real, nₜ::Real, N::Integer, h::Real, method::String="RK4")

Resuelve dos ecuaciones diferenciales que representan una reacción de dimerización.


# Argumentos
- `t₀::Real`: Tiempo inicial.
- `χ₀::Vector{Float64}`: Condiciones iniciales para las especies químicas.
- `k₁::Real`, `k₂::Real`: Constantes de la reacción.
- `nₜ::Real`: Concentración total de las especies.
- `N::Integer`: Número de pasos en el método numérico.
- `h::Real`: Tamaño del paso en el método numérico.
- `method::String="RK4"`: Método numérico para resolver las ecuaciones diferenciales. Puede ser "Euler", "MidPoint", "Ralston" o "RK4".

# Salida
- `χ₁`, `χ₂`: Soluciones de las ecuaciones diferenciales para cada especie.
- `t`: Vector de tiempos.

# Ejemplos
```julia
χ₁, χ₂, t = Dimerization(0.0, [1.0, 0.0], 1.0, 1.0, 1.0, 100, 0.01, "RK4")
"""


function dimerization(t₀::Real, χ₀₁::Real, χ₀₂::Real, k₁::Real, k₂::Real, nₜ::Real=1.0, N::Integer=10000, h::Real=0.01, method::String="RK4")
    # Definimos las funciones para las derivadas
    f₁ = (t, χ) -> -2*k₁*nₜ*(χ^2) + k₂*(1 - χ)
    f₂ = (t, χ) -> 2*k₁*nₜ*((1 - χ)^2) - k₂*χ

    # Resolvemos las ecuaciones diferenciales para cada especie
    χ₁, χ₂, t₁ = sol_two_ODE(f₁, f₂, χ₀₁, χ₀₂, t₀, h, N, method=method)
    return χ₁, χ₂, t₁
end

"""
    dimerization_stationary(α₁::Real, α₂::Real, t::Vector{T}, nₜ::Real=1.0) where {T<:Real}

Calcula las soluciones estacionarias para un sistema de dimerización.

# Argumentos
- `α₁::Real`: Parámetro de la ecuación.
- `α₂::Real`: Parámetro de la ecuación.
- `t::Vector{T}`: Vector de tiempo.
- `nₜ::Real=1.0`: Parámetro de normalización.

# Devoluciones
- `valid_solutions::Vector{T}`: Vector con las dos soluciones válidas en el rango [0,1].
"""
function dimerization_stationary(α₁::Real, α₂::Real, t::Vector{T}, nₜ::Real=1.0) where {T<:Real}
    Γ = α₂ / (2 * α₁ * nₜ)
    sqrt_term = sqrt(Γ^2 + 4 * Γ)

    χ₁ˢˢ_vals = [(-Γ + sqrt_term) / 2, (-Γ - sqrt_term) / 2]
    χ₂ˢˢ_vals = [((2 + Γ) + sqrt_term) / 2, ((2 + Γ) - sqrt_term) / 2]

    χ₁₁ˢˢ, χ₁₂ˢˢ = fill.(χ₁ˢˢ_vals, length(t))
    χ₂₁ˢˢ, χ₂₂ˢˢ = fill.(χ₂ˢˢ_vals, length(t))

    solutions = [χ₁₁ˢˢ, χ₁₂ˢˢ, χ₂₁ˢˢ, χ₂₂ˢˢ]
    valid_solutions=[]
    for solution in solutions
        if check_values_in_range(solution)
            push!(valid_solutions, solution)
        end
    end
#    if length(valid_solutions) == 2
#        println("Las dos soluciones en el rango [0,1] son: ", valid_solutions)
    return valid_solutions

end


