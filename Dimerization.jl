#About: This code search solve the two first moments for the distribution probability of a stochastic process defined by the NonEquilibrium theory of Onsager-Machlup to describe the dimerization, and solutions for the Keizer theory (canonical molecular process)
#Author: Jan Carlo Alvarez Centeno & Ricardo Peredo Ortiz
# Institution: Intituto de Física UASLP
#------------Packages---------------
using Plots
using RecipesBase
using DelimitedFiles

#-----------------------------------
#--------------------------------Method to integrate a differential equation --------------------------------
"""
Here we define different methods to numerical evolution
inputs:
------
	f :: Function # Function of the variables t and x to be evalauted
	x₀ :: Real # Initial x value
	t₀ :: Real # Initial t value
	h :: Real # differential dt element
	N :: Integer # Number of t-steps
	method = "RK4" # Method used to evaluate the ODE
"""
function sol_diff_eq_ord(f :: Function, x₀ :: Real, t₀ :: Real, h :: Real, N :: Integer; method = "RK4")
	t = collect(t₀:h:t₀ + N*h)
	x = zeros(length(t))
	x[1] = x₀
	if method == "Euler"
		for n in 2:length(t)
			x[n] = x[n-1] + h*f(t[n-1], x[n-1])
		end
		return x, t
	elseif method == "MidPoint"
		for n in 2:length(t)
			x[n] = x[n-1] + h*f(t[n-1] + 0.5*h, x[n-1] + 0.5*h*f(t[n-1], x[n-1]))
		end
		return x, t
	elseif method == "Ralston"
		for n in 2:length(t)
			k₁ = f(t[n-1], x[n-1])
			k₂ = f(t[n-1] + (2/3)*h, x[n-1] + (2/3)*h*k₁)
			x[n] = x[n-1] + h*(0.25*k₁ + 0.75*k₂)
		end
		return x, t
		return x, t
	elseif method == "RK4"
		for n in 2:length(t)
			k₁ = f(t[n-1], x[n-1])
			k₂ = f(t[n-1] + 0.5*h, x[n-1] + 0.5*h*k₁)
			k₃ = f(t[n-1] + 0.5*h, x[n-1] + 0.5*h*k₂)
			k₄ = f(t[n-1] + h, x[n-1] + h*k₃)
			x[n] = x[n-1] + h*(k₁ + 2*k₂ + 2*k₃ + k₄)/6
		end
		return x, t
	end
end;
#--------------------------------------------------------------------------------------------------------------
#--------------Functions to integrate---------------------
# The varible is the molar fraction χ₁+χ₂=1.0, where χ₁=n₁/nₜ with n₁ is the volumentric concentration [1/V] and nₜ is the total concentration. So χ is adimentional.
function Dimerization_n1(t, χ₁, k₁ ,k₂ ,nₜ)
    return -2*k₁*nₜ*(χ₁^2)+k₂-k₂*χ₁
end
function Dimerization_n2(t, χ₂, k₁ ,k₂ ,nₜ)
    return 2*k₁*nₜ*(1.0-χ₂)^2-k₂*χ₂
end
#---------------------------------------------------------
#------------Parameters------------
χ₁₀ = 1.0 # Initial condition for concentration
χ₂₀ = 0.0
nₜ= 1.0 # Total volumetric concentration
t₀ = 0.0 # Initial time
h = 0.01 # Length for the step time
N_steps = 10000 # numbre of steps in the time 
α₁ = 0.9 # Rate constant forward
α₂ = 0.1 # Rate constant reverse
λ = α₁+α₂
#----------------------------------
#-------------Redefine the functions-------------
Dim₁(t,χ₁)=Dimerization_n1(t,χ₁,α₁,α₂,nₜ)
Dim₂(t,χ₂)=Dimerization_n2(t,χ₂,α₁,α₂,nₜ)
#------------------------------------------------
#-----------------------------------For_Plot--------------------------------------- 
plot(title = "Dimerización k₁="*string(α₁)*",k₂="*string(α₂), xlabel = "t-Tiempo", ylabel = "n(t) - Concentración [1/V]", legend = :topright)
χ₁, t₁ = sol_diff_eq_ord(Dim₁, χ₁₀, t₀, h, N_steps)
plot!(t₁, χ₁, label = "χ₁", color="darkgreen", linewidth=3)
χ₂, t₂ = sol_diff_eq_ord(Dim₂, χ₂₀, t₀, h, N_steps)
plot!(t₂, χ₂, label = "χ₂", xlims = (0,6), color="violet", linewidth=3)
Γ=α₂/(2*α₁*nₜ)
χ₁ˢˢ_val_1=(-Γ+sqrt(Γ^2+4*Γ))/2
χ₁ˢˢ_val_2=(-Γ-sqrt(Γ^2+4*Γ))/2
χ₂ˢˢ_val_1=((2+Γ)+sqrt(Γ^2+4*Γ))/2
χ₂ˢˢ_val_2=((2+Γ)-sqrt(Γ^2+4*Γ))/2
χ₁ˢˢ=fill(χ₁ˢˢ_val_1,length(t₁))
χ₂ˢˢ=fill(χ₂ˢˢ_val_2,length(t₂))
plot!(t₁,χ₁ˢˢ, linewidth=2, label="n₁ - estacionaria ", color="green", line=:dash)
plot!(t₂,χ₂ˢˢ, linewidth=2, label="n₂ - estacionaria ", color="violetred", line=:dash)
savefig("Dimerizacion2.png")
#----------------------------------------------------------------------------------- 
#----------To save the data-------------
open("Dimerzacion.csv","w") do file
    writedlm(file, [χ₁ t₁ χ₂ t₂], ';')
end
#---------------------------------------
function sol_diff_eq_ord_1par(f :: Function, c::Array, x₀ :: Real, t₀ :: Real, h :: Real, N :: Integer; method = "RK4")
	t = collect(t₀:h:t₀ + N*h)
	x = zeros(length(t))
	x[1] = x₀
	if method == "Euler"
		for n in 2:length(t)
			x[n] = x[n-1] + h*f(t[n-1], x[n-1],c[n-1])
		end
		return x, t
	elseif method == "MidPoint"
		for n in 2:length(t)
			x[n] = x[n-1] + h*f(t[n-1] + 0.5*h, x[n-1] + 0.5*h*f(t[n-1], x[n-1],c[n-1]),c[n-1])
		end
		return x, t
	elseif method == "Ralston"
		for n in 2:length(t)
			k₁ = f(t[n-1], x[n-1],c[n-1])
			k₂ = f(t[n-1] + (2/3)*h, x[n-1] + (2/3)*h*k₁,c[n-1])
			x[n] = x[n-1] + h*(0.25*k₁ + 0.75*k₂)
		end
		return x, t
		return x, t
	elseif method == "RK4"
		for n in 2:length(t)
			k₁ = f(t[n-1], x[n-1],c[n-1])
			k₂ = f(t[n-1] + 0.5*h, x[n-1] + 0.5*h*k₁,c[n-1])
			k₃ = f(t[n-1] + 0.5*h, x[n-1] + 0.5*h*k₂,c[n-1])
			k₄ = f(t[n-1] + h, x[n-1] + h*k₃,c[n-1])
			x[n] = x[n-1] + h*(k₁ + 2*k₂ + 2*k₃ + k₄)/6
		end
		return x, t
	end
end; 
#-----------------Variance Equations------------------
σˢˢ=[(2*χ₁[i]*(1.0-χ₁[i])/(2-χ₁[i])) for i in 1:length(t₁)]
H=[-(4*α₁*χ₁[i]*nₜ+α₂) for i in 1:length(t₁)]
function variance_dim_OMNE(t,σ,χ;k_1=α₁,k_2=α₂,nₜ=nₜ)
	H=-(4*k_1*χ*nₜ+k_2)
	σˢˢ=(2*χ*(1.0-χ))/(2-χ)
	return 2*H*(σ-σˢˢ)
end
#function variance_iso_Keizer()
#end
#-----------------------------------------------------
#------------------Obtain the variance---------------------
σ₀=0.0
σ_dim, t_dim=sol_diff_eq_ord_1par(variance_dim_OMNE, χ₁, σ₀, t₀, h, N_steps)
σˢˢ₁_val=(2*χ₁ˢˢ_val_1*(1-χ₁ˢˢ_val_1))/(2-χ₁ˢˢ_val_1)
σˢˢ₁=fill(σˢˢ₁_val,length(t_dim))
#σˢˢ₂=fill(,length(tᵢₛₒ))
#----------------------------------------------------------
function Trapezoidal_rule(f::Function,x_min::Real,x_max::Real,N_points::Integer)
    Δx=(x_max-x_min)/N_points #This value is the length of the steps for integration
    sum = 0.0 #Initialice value for the sum
    for k=1:N_points-1
		step_1=n_min+(k-1)*Δx
		step_2=n_min+(k)*Δx
        f1 = f(step_1)
        f2 = f(step_2)
        sum = sum + Δx*(f1+f2)/2
    end
    return sum
end
#--------------------------------------------------------------------------------------------
#--------------------------- Trapezoidal rule with one array parameter ----------------------------
function Trapezoidal_rule_1_parameter(f::Function,c::Array,t::Array,h::Real)
	# In the case of Isomerization c is the concentration and t is the time value
	resultado=fill(0.0,length(t))
	for i in 2:length(t)
		sum=0.0
		Δx=h
		for j in 2:i
			f1 = f(c[j-1])
			f2 = f(c[j])
			sum = sum + Δx*(f1+f2)/2
		end	
		resultado[i]=sum
	end
	return resultado
end
function Trapezoidal_rule_2_parameter(f::Function,c::Array,μ::Array,t::Array,h::Real)
	# In the case of Isomerization c is the concentration and t is the time value
	resultado=fill(0.0,length(t))
	for i in 2:length(t)
		sum=0.0
		Δx=h
		for j in 2:i
			f1 = μ[j-1]*f(c[j-1])
			f2 = μ[j]*f(c[j])
			sum = sum + Δx*(f1+f2)/2
		end	
		resultado[i]=sum
	end
	return resultado
end
#--------------------------------------------------------------------------------------------------------
#--------------------Gamma and H function for Keizer Theory------------------------
function gamma_function_dim(χ₁::Real,k₁::Real,k₂::Real,nₜ::Real)
	return 4*(k₁*(χ₁^2)+0.5*k₂*(1.0-χ₁))	
end
function H_fluctuation_dim(χ₁::Real, k₁::Real, k₂::Real,nₜ::Real)
	return -2*(4*k₁*χ₁*nₜ+k₂)
end
γₖ(x_1)=gamma_function_dim(x_1,α₁,α₂,nₜ)
Hₖ(x_1)=H_fluctuation_dim(x_1,α₁,α₂,nₜ)
Int_do=Trapezoidal_rule_1_parameter(Hₖ,χ₁,t₁,h)
#----------------------------------------------------------------------------------
#-----------------------Variance Keizer--------------------------------------------
μ₁=[exp(Int_do[i]) for i in 1:length(Int_do)]
μ₂=[exp(-Int_do[i]) for i in 1:length(Int_do)]
Int_σ=Trapezoidal_rule_2_parameter(γₖ,χ₁,μ₂,t₁,h)
σₖ=[μ₁[i]*Int_σ[i] for i in 1:length(t₁)]
#----------------------------------------------------------------------------------
#----------------Data-----------------
#open("Variance_Dim.csv","w") do file
#    writedlm(file, [σ_dim σˢˢ H t_dim], ';')
#end
#-------------------------------------
#------------------Plot----------------
plot(t_dim,σ_dim, xlims=(0.0,10.0), title="Dimerización k₁="*string(α₁)*" k₂="*string(α₂), linewidth=3, label="OM-NE")
plot!(t_dim, σˢˢ₁, label="σˢˢ", linewidth=2, line=:dash)
plot!(t_dim, σˢˢ, label="σˢˢ(t)")
plot!(t_dim, σₖ, label="Keizer",xlabel="Tiempo-t",ylabel="Varianza-σ(t)")
savefig("Variance_Dim.png")
#--------------------------------------