using Plots

"""
    plot_dimerization(t, χ₁, χ₂, χ₁ˢˢ, χ₂ˢˢ)

Función para graficar las soluciones de la dimerización y las soluciones estacionarias.

# Argumentos
- `t`: Vector de tiempo.
- `χ₁`: Vector de soluciones de la dimerización χ₁.
- `χ₂`: Vector de soluciones de la dimerización χ₂.
- `χ₁ˢˢ`: Vector de soluciones estacionarias χ₁ˢˢ.
- `χ₂ˢˢ`: Vector de soluciones estacionarias χ₂ˢˢ.

# Ejemplo
```julia
t = 0:0.1:10
χ₁ = sin.(t)
χ₂ = cos.(t)
χ₁ˢˢ = sin.(t) .+ 0.5
χ₂ˢˢ = cos.(t) .+ 0.5
plot_dimerization(t, χ₁, χ₂, χ₁ˢˢ, χ₂ˢˢ)

"""

function plot_dimerization(t, χ₁, χ₂, χ₁ˢˢ, χ₂ˢˢ)
    # Crear una nueva figura
    plt = plot()

    # Agregar las soluciones de la dimerización
    plot!(plt, t, χ₁, label="χ₁(t)", linewidth=2)
    plot!(plt, t, χ₂, label="χ₂(t)", linewidth=2)

    # Agregar las soluciones estacionarias
    plot!(plt, t, χ₁ˢˢ, label="χ₁ˢˢ", linewidth=2, linestyle=:dash)
    plot!(plt, t, χ₂ˢˢ, label="χ₂ˢˢ", linewidth=2, linestyle=:dash)

    # Configurar los límites del eje x
    xlims!(plt, 0, 10)

    # Configurar los títulos de los ejes y la gráfica
    title!(plt, "Soluciones de la Dimerización y Soluciones Estacionarias")
    xlabel!(plt, "Tiempo")
    ylabel!(plt, "Soluciones")

    # Mostrar la gráfica
    display(plt)
    savefig("images/Dimerizacion_27_05_24.png")
end

