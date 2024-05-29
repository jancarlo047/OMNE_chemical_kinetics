#include("generalfunctions.jl")
include("generalChemFunc.jl")
using Plots


χ₁₀ = 1.0 # Initial condition for concentration
χ₂₀ = 0.0
t₀ = 0.0 # Initial time
α₁ = 0.9 # Rate constant forward
α₂ = 0.1 # Rate constant reverse
λ = α₁+α₂
χ₁,χ₂,t= dimerization(t₀,χ₁₀,χ₂₀,α₁,α₂)
χ₁ˢˢ, χ₂ˢˢ = dimerization_stationary(α₁, α₂, collect(t))

plot_dimerization(t, χ₁, χ₂, χ₁ˢˢ, χ₂ˢˢ)
