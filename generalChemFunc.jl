include("generalfunctions.jl")

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

function Dimerization(t₀::Real, χ₀::Real, k₁::Real, k₂::Real, nₜ::Real=1.0, N::Integer=10000, h::Real=0.01, method::String="RK4")
    # Definimos las funciones para las derivadas
    f₁ = (t, χ) -> -2*k₁*nₜ*(χ^2) + k₂*(1 - χ)
    f₂ = (t, χ) -> 2*k₁*nₜ*((1 - χ)^2) - k₂*χ

    # Resolvemos las ecuaciones diferenciales para cada especie
    χ₁, t₁ = sol_diff_eq_ord(f₁, χ₀, t₀, h, N, method=method)
    χ₂, t₂ = sol_diff_eq_ord(f₂, 1.0-χ₀, t₀, h, N, method=method)
#    print(string(typeof(χ₁)) * "\n")
#    print(string(typeof(χ₂)) * "\n")
#    print(string(typeof(t₁)) * "\n")
    return χ₁, χ₂, t₁
end