using Distributions
using Combinatorics
using LinearAlgebra
using Statistics
using Random

function calc_shapley(Z; vals = 1:size(Z)[2]-1, with_cov = false)
  d, n = size(Z)[2]-1, size(Z)[1]

  # Equation (2) (taking care of 0-indexing)
  Cₙ_ = cor(Z)
  Cₙ(u) = [Cₙ_[i,j] for i in u .+ 1, j in u .+ 1]

  # Equation (4) (taking care of empty s here)
  R²(s) = (length(s) > 0) ? 1 - det(Cₙ(vcat(0,s)))/det(Cₙ(s)) : 0

  # Equation (9) (pre-compute ω_ for efficiency)
  ω_ = [factorial(i)/factorial(d, d-i-1) for i in 0:(d-1)]
  ω(s) = ω_[length(s) + 1]
  S(j) = deleteat!(collect(1:d), j)
  V(j) = sum([ω(s)*(R²(vcat(j,s)) - R²(s)) for s in powerset(S(j))])

  # Remark 4
  μZ = [mean(Z[:,j+1]) for j in 0:d]
  invCov = inv(cov(Z))
  mahalanobis_sq = [(Z[i,:] .- μZ)'*invCov*(Z[i,:] .- μZ) for i in 1:n]
  κ = sum( mahalanobis_sq.^2 )/(n*(d+1)*(d+3))

  # Equation (6) (pre-compute acovζ_ for efficiency)
  function acovζ(g,h,j,k)::Float64
    g += 1; h += 1; j += 1; k += 1 # taking care of 0-indexing
    quad_part = Cₙ_[g,j]^2 + Cₙ_[h,j]^2 + Cₙ_[g,k]^2 + Cₙ_[h,k]^2
    κ*(Cₙ_[g,h]*Cₙ_[j,k]*(quad_part)/2
      + Cₙ_[g,j]*Cₙ_[h,k] + Cₙ_[g,k]*Cₙ_[h,j]
      - Cₙ_[g,h]*(Cₙ_[h,j]*Cₙ_[h,k] + Cₙ_[g,j]*Cₙ_[g,k])
      - Cₙ_[j,k]*(Cₙ_[g,j]*Cₙ_[h,j] + Cₙ_[g,k]*Cₙ_[h,k]))
  end
  acovζ_ = [acovζ(g,h,j,k) for g in 0:d, h in 0:d, j in 0:d, k in 0:d]

  # Equation (7) (taking care of 0-indexing in acovδ)
  function R(u)
    if (length(u) == 0) return 0 end
    det(Cₙ(u))*inv(Cₙ(u))
  end
  function acovδ(u, v)
    if (length(u) == 0 || length(v) == 0) return 0 end
    rᵤ, rᵥ = R(u), R(v)
    lu,lv = 1:length(u), 1:length(v)
    rᵤrᵥacovζ(g,h,j,k)=rᵤ[g,h]*rᵥ[j,k]*acovζ_[u[g].+1,u[h].+1,v[j].+1,v[k].+1]
    sum([rᵤrᵥacovζ(g,h,j,k) for g in lu, h in lu, j in lv, k in lv])
  end

  # Equation (8)
  function acovλ(s,t)::Float64
    if (length(s) == 0 || length(t) == 0) return 0 end
    if ( s == t ) return 4 .* R²(s) .* (1 .- R²(s)) .^2 end
    CₙS, CₙT, CₙS0, CₙT0 = Cₙ(s), Cₙ(t), Cₙ(vcat(0,s)), Cₙ(vcat(0,t))
    detCₙS, detCₙT, detCₙS0, detCₙT0 = det(CₙS), det(CₙT), det(CₙS0), det(CₙT0)
    (1/(detCₙS*detCₙT)*acovδ(vcat(0,s),vcat(0,t))
      + detCₙS0*detCₙT0/(detCₙS^2*detCₙT^2)*acovδ(s,t)
      - detCₙS0/(detCₙS^2*detCₙT)*acovδ(s,vcat(0,t))
      - detCₙT0/(detCₙS*detCₙT^2)*acovδ(vcat(0,s),t))
  end

  # Theorem 1
  function acovξ(j,k)::Float64
    powSⱼ, powSₖ  = powerset(S(j)), powerset(S(k))
    sum([ω(s)*ω(t)*(
         acovλ(vcat(j,s), vcat(k,t)) +
         acovλ(s,t) -
         acovλ(s, vcat(k,t)) -
         acovλ(vcat(j,s), t)) for s in powSⱼ, t in powSₖ])
  end

  # Calculate all the shapley values using Equation (9)
  shapley = map( x -> V(x), vals )

  if (!with_cov)
    # Calculate shapley variances using Theorem 1
    shapley_var = map( j -> acovξ(j,j), vals )
    return shapley, shapley_var
  else
    # Calculate covariances using Theorem 1
    shapley_cov = [[acovξ(i,j) for i in filter(i -> i <= j, vals)] for j in vals]
    return shapley, shapley_cov
  end
end

# Calculate confidence intervals using Theorem 1
function calc_CIs(shap, avar, n; α = 0.05)
  Φ⁻¹ = quantile.(Normal(), 1-α/2)
  Δ = Φ⁻¹ .* map(sqrt, map(abs, avar) ./ n)
  CIℓ = shap .- Δ
  CIu = shap .+ Δ
  return CIℓ, CIu
end


# Since we know the exact correlation matrix, we can calculate
# the true shapley values vⱼ = [ limit of E(Vⱼ(Zₙ)) as n → ∞ ]
function true_shapley(dist; vals = 1)
  # The exact correlation matrix of the simulated distribution
  R_ = cor(dist)

  # Using Equation (4) and Slutsky's Theorem, as n → ∞,
  # the expected value of R² is asymptotically unbiased.
  R(u) = [R_[i,j] for i in u .+ 1, j in u .+ 1]
  ER²(s) = (length(s) > 0) ? 1 - det(R(vcat(0,s)))/det(R(s)) : 0

  # Expected value of Equation (9), as n → ∞
  # (pre-compute ω_ for efficiency)
  ω_ = [factorial(i)/factorial(d, d-i-1) for i in 0:(d-1)]
  ω(s) = ω_[length(s) + 1]
  S(j) = deleteat!(collect(1:d), j)
  v(j) = sum([ω(s)*(ER²(vcat(j,s)) - ER²(s)) for s in powerset(S(j))])
  map( x -> v(x), vals )
end

# This function, used for bootstrap,
# just calculates the shapley values
# and does not bother with the variances,
# so it's just the first half of calc_shapley
function just_shapley(Z; vals = 1:size(Z)[2]-1)
  d, n = size(Z)[2]-1, size(Z)[1]

  # Equation (2) (taking care of 0-indexing)
  Cₙ_ = cor(Z)
  Cₙ(u) = [Cₙ_[i,j] for i in u .+ 1, j in u .+ 1]

  # Equation (4) (taking care of empty s here)
  R²(s) = (length(s) > 0) ? 1 - det(Cₙ(vcat(0,s)))/det(Cₙ(s)) : 0

  # Equation (9) (pre-compute ω_ for efficiency)
  ω_ = [factorial(i)/factorial(d, d-i-1) for i in 0:(d-1)]
  ω(s) = ω_[length(s) + 1]
  S(j) = deleteat!(collect(1:d), j)
  V(j) = sum([ω(s)*(R²(vcat(j,s)) - R²(s)) for s in powerset(S(j))])

  # Calculate all the shapley values using Equation (9)
  return map( x -> V(x), vals )
end













# #### TESTING
# # Parameters
# Random.seed!(0)
# d, n = 2, 1000
# X = randn( n, d )
# y = X * collect(0:2:(2d-2)) + randn( n )
# Z = hcat(y,X)
# @time s, v = calc_shapley(Z)
# @time s, v = calc_shapley(Z, with_cov = true)
#
# # Check that compiler figures out all the types
# @code_warntype calc_shapley(Z)
# # Benchmarking
# using BenchmarkTools
# @benchmark calc_shapley(Z)
