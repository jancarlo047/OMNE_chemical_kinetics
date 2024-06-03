"""
    check_values_in_range(vector::Vector{Float64})

Revisa si todos los valores en un vector están en el rango [0,1].

# Argumentos
- `vector::Vector{Float64}`: Un vector de números de punto flotante.

# Devoluciones
- `Bool`: Devuelve `true` si todos los valores están en el rango [0,1], de lo contrario devuelve `false`.

# Ejemplos
```julia
julia> check_values_in_range([0.1, 0.2, 0.3])
true

julia> check_values_in_range([0.1, 1.2, 0.3])
false"""


function check_values_in_range(vector::Vector{Float64})
    return all(x -> 0 <= x <= 1, vector)
end
