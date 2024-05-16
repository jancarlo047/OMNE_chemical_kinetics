#About: This code search solve the two first moments for the distribution probability of a stochastic process defined by the Non-Equilibrium Theory of Fluctuations by Onsager-Machlup process to describe the Isomerization reaction
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
function Isomerization_n1(t, n₁, k₁ ,k₂, nₜ)
    return -(k₁+k₂)*n₁+k₂*nₜ
end
function Isomerization_n2(t, n₂, k₁ ,k₂, nₜ)
    return -(k₁+k₂)*n₂+k₁*nₜ
end
#---------------------------------------------------------
#------------Parameters------------
n₁₀ = 1.0 # Initial condition for concentration
nₜ = n₁₀
t₀ = 0.0 # Initial time
h = 0.01 # Length for the step time
N_steps = 10000 # numbre of steps in the time 
α₁ = 0.1 # Rate constant forward
α₂ = 0.9 # Rate constant reverse
λ = α₁+α₂
#----------------------------------
#-------------Redefine the functions-------------
Iso₁(t,n₁)=Isomerization_n1(t,n₁,α₁,α₂,nₜ)
Iso₂(t,n₂)=Isomerization_n2(t,n₂,α₁,α₂,nₜ)
#------------------------------------------------
#-----------------------------------For_Plot--------------------------------------- 
plot(title = "Isomerización k₁="*string(α₁)*",k₂="*string(α₂), xlabel = "t-Tiempo", ylabel = "n(t) - Concentración [1/V]", legend = :topright)
n₁, t₁ = sol_diff_eq_ord(Iso₁, n₁₀, t₀, h, N_steps)
plot!(t₁, n₁, label = "n₁", color="darkgreen", linewidth=3)
n₂, t₂ = sol_diff_eq_ord(Iso₂, nₜ-n₁₀, t₀, h, N_steps)
plot!(t₂, n₂, label = "n₂", xlims = (0,6), color="violet", linewidth=3)
ss_n1=fill(α₂*nₜ/λ,length(t₁))
ss_n2=fill(α₁*nₜ/λ,length(t₂))
plot!(t₁,ss_n1, linewidth=2, label="n₁ - estacionaria ", color="green", line=:dash)
plot!(t₂,ss_n2, linewidth=2, label="n₂ - estacionaria ", color="violetred", line=:dash)
savefig("Isomerizacion.png")
#----------------------------------------------------------------------------------- 
#-----------------------Obtain data for save-------------------------
#
#----------To save the data-------------
open("Isomerizacion.csv","w") do file
    writedlm(file, [n₁ t₁ n₂ t₂], ';')
end
#---------------------------------------
#----------------Numerical Solution of an differential equation with parameter----------------
function sol_diff_eq_ord_1_parameter(f :: Function, c::Array, x₀ :: Real, t₀ :: Real, h :: Real, N :: Integer; method = "RK4")
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
#------------------------------------------------------------------------------------- 
#-----------------Variance Equations------------------
σˢˢₜ=[n₁[i]*n₂[i]/nₜ for i in 1:length(t₁)]
function variance_iso_OMNE(t,σ,n;k_1=α₁,k_2=α₂,nₜ=nₜ)
	λ=k_1+k_2
	return -2*λ*(σ-n*(1.0-n))
end
function variance_Keizer(σ₀::Real, k₁::Real, k₂::Real,  aₜ::Real, t::Real)
	λ=k₁+k₂
	return (k₁*aₜ/(λ*λ))*((k₁-k₂)*(exp(-λ*t)-exp(-2*λ*t))+k₂*(1-exp(-2*λ*t)))
end
#-----------------------------------------------------
#------------------Obtain the variance---------------------
σ₀=0.0
σᵢₛₒ, tᵢₛₒ=sol_diff_eq_ord_1_parameter(variance_iso_OMNE, n₁, σ₀, t₀, h, N_steps)
σˢˢ=fill((n₁₀*α₁*α₂)/(λ^2),length(tᵢₛₒ))
f1(t)=variance_Keizer(0.0, α₁, α₂, nₜ,t)
#----------------------------------------------------------
#------------------------------------ Trapezoidal rule --------------------------------------
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
	resultado=fill(0.0,length(tᵢₛₒ))
	for i in 2:length(tᵢₛₒ)
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
	resultado=fill(0.0,length(tᵢₛₒ))
	for i in 2:length(tᵢₛₒ)
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
function gamma_function_iso(n₁,k₁::Real,k₂::Real,nₜ::Real)
	return k₁*n₁+k₂*(nₜ-n₁)	
end
function H_fluctuation_iso(n₁::Real, k₁::Real, k₂::Real)
	return -2(k₁+k₂)
end
γₖ(n_1)=gamma_function_iso(n_1,α₁,α₂,nₜ)
Hₖ(n_1)=H_fluctuation_iso(n_1,α₁,α₂)
Int_do=Trapezoidal_rule_1_parameter(Hₖ,n₁,t₁,h)
#----------------------------------------------------------------------------------
#-----------------------Variance Keizer--------------------------------------------
μ₁=[exp(Int_do[i]) for i in 1:length(Int_do)]
μ₂=[exp(-Int_do[i]) for i in 1:length(Int_do)]
Int_σ=Trapezoidal_rule_2_parameter(γₖ,n₁,μ₂,t₁,h)
σₖ=[μ₁[i]*Int_σ[i] for i in 1:length(t₁)]
#----------------------------------------------------------------------------------
#----------------Data-----------------
#open("Variance_Iso.csv","w") do file
#    writedlm(file, [σᵢₛₒ tᵢₛₒ], ';')
#end
#-------------------------------------
#------------------Plot----------------
plot(tᵢₛₒ,σᵢₛₒ, xlims=(0.0,10.0), title="Isomerización k₁="*string(α₁)*" k₂="*string(α₂), linewidth=3, label="OM-NE")
plot!(tᵢₛₒ, σˢˢ, label="σˢˢ=x₀k₁k₂/λ²", linewidth=2, line=:dash)
plot!(tᵢₛₒ, σˢˢₜ, label="Estacionaria")
#plot!(f1,0,50, label="Keizer & McQuarrie", linewidth=3)
plot!(t₁,σₖ,linewidth=3, label="Keizer",xlabel="Tiempo-t",ylabel="Varianza-σ(t)")
savefig("Variance_Iso.png")
#--------------------------------------
