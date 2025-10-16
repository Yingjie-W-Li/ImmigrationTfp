"""
Parameters of the auxilliary model for TFP. I call this an auxilliary model because
it could be repeated on model-simulated data.
"""
struct AuxParameters{T1 <: Real}

    ρ::T1                                       # CES parameter
    θ::T1                                       # Capital share
    αᶠ::T1                                      # Absolute advantage foreign
    αᵈ::T1                                      # Absolute advantage domestic
    Inter::T1                                   # Intercept
    δ::Vector{T1}                               # State fixed effects
    ξ::Vector{T1}                               # Time fixed effects
    ζᶠ::T1                                      # Comp advantage - foreign
    ζᵈ::T1                                      # Comp advantage - domestic

end

"""
Constructor for the AuxParameters type.
    - df should be sorted, state then year. Thus, reference levels for fixed effects are Alabama in 1994
"""
function AuxParameters(;
    ρ::T1 = 1.,
    θ::T1 = 0.3,
    αᶠ::T1 = 1.,
    αᵈ::T1 = 1.,
    Inter::T1 = 1.,
    ζᶠ::T1 = 1.,
    ζᵈ::T1 = 1.,
    df::DataFrame = Wide,
    N::Int64 = length(unique(df[:,:statefip])),             # Number of years
    T::Int64 = length(unique(df[:,:year])),                 # Number of units
    δ::Vector{T1} = vcat(0., ones(N - 1)),
    ξ::Vector{T1} = vcat(0., ones(T - 1)),
    ) where{T1 <: Real}                
    
    return AuxParameters{T1}(ρ, θ, αᶠ, αᵈ, Inter, δ, ξ, ζᶠ, ζᵈ)

end

"""
Calculate the cutoff task.
"""
function 𝒯(p::AuxParameters; df::DataFrame = Wide)

    (; ζᶠ, ζᵈ, ρ, αᶠ, αᵈ) = p
    
    # Define parameters for readability
    b = ρ/(1 - ρ)
    α = αᶠ/αᵈ
    ζ = (1 + ζᵈ * b)/(1 + ζᶠ * b)
    wᵈ = df[:, :Wage00]
    wᶠ = df[:, :Wage01]


    return ((wᵈ ./ wᶠ) * α * ζ^(-(1/b))).^(1/(ζᵈ - ζᶠ))

end

"""
The residual function. 
    - df should be sorted state and then year
    - The first element of λ and δ should be 0 to avoid collinearity
"""
function Residual(p::AuxParameters; df::DataFrame = Wide)
    
    (; δ, ξ, ζᶠ, ζᵈ, θ, ρ, αᶠ, αᵈ, Inter) = p

    T = length(unique(df[:,:year]))                        # Number of years
    N = length(unique(df[:,:statefip]))                    # Number of units
    TFE = repeat(Matrix{Float64}(I, T, T), N)              # Time FE mat

    # Set up the state fixed effects matrix
    SFE = Matrix{Float64}(undef, 0, N)
    for c in 1:N
        SFE = vcat(SFE, [j == c ? 1. : 0. for i in 1:T, j in 1:N])
    end

    # Calculate task shares
    b = ρ/(1 - ρ)
    T_cal = max.(min.(𝒯(p; df = df), 1.),0.)
    λ = T_cal.^(1 + ζᶠ* b)/(1 .+ T_cal.^(1 + ζᶠ * b) .- T_cal.^(1. + ζᵈ * b))

    # Calculate each part of the production function
    Z = (1 .+ T_cal.^(1 + ζᶠ * b) .- T_cal.^(1 + ζᵈ* b)).^(1/b)
    K = df[:,:K]
    F = df[:,:BodiesSupplied01]
    D = df[:, :BodiesSupplied00]
    L = (λ.^(1 - ρ) * (αᶠ * F).^ρ + (1 .- λ).^(1 - ρ) * (αᵈ * D).^ρ).^(1/ρ)
    Y = df[:, :Y]

    return log.(Y) - (Inter .+ SFE * δ + TFE * ξ + θ * log.(K) + (1 - θ) * log.(Z .* L)), Z, L
    
end

"""
The sum of squared errors for state level production function.
"""
function SSE(x::Vector{T1}; df::DataFrame = Wide) where{T1 <: Real}

    N = length(unique(df[:, :statefip]))
    T = length(unique(df[:, :year]))
    ρ, θ, αᶠ, αᵈ, Inter, δ, ξ, ζᶠ, ζᵈ = x[1], x[2], x[3], x[4], x[5], vcat(0., x[6: 5 + N - 1]), vcat(0., x[5 + N : 5 + N - 1 + T - 1]), x[end - 1], x[end]

    p = AuxParameters(ρ, θ, αᶠ, αᵈ, Inter, δ, ξ, ζᶠ, ζᵈ)

    vals = Residual(p; df = df)

    return norm(vals[1]).^2

end


