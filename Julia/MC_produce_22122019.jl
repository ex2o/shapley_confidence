using Distributions
using LinearAlgebra
using Statistics
using Random
using BenchmarkTools
using CSV, DataFrames
using HypothesisTests
include("Shapley.jl")

### HELPER FUNCTIONS

# (d+1)x(d+1) correlation matrix with uniform correlations c
M0(c,d) = [Float64(c) + (i == j)*(1-Float64(c)) for i in 1:(d+1), j in 1:(d+1)]
wish(c, d, df) = Wishart(df, M0(c,d))
W0(wdist) = rand(wdist, 1)[1]


# Get summary statistics about CI width and coverage
function CI_measures(within, width)
  if (mean(within) == 0)
    covrCI = (0,0)
  elseif (mean(within) == 1)
    covrCI = (1,1)
  else
    covrCI = confint(BinomialTest(within))
  end
  widthCI = quantile!(width, [0.025, 0.975])
  return mean(within), covrCI[1], covrCI[2], mean(width),
         std(width), widthCI[1], widthCI[2]
end

# Function to save results to csv
function save_sim_csv(sim, type, d_)
  CSV.write("results/$d_/sim$type/coverage.csv", DataFrame(sim[1]), writeheader=false)
  CSV.write("results/$d_/sim$type/covrgCI_L.csv", DataFrame(sim[2]), writeheader=false)
  CSV.write("results/$d_/sim$type/covrgCI_U.csv", DataFrame(sim[3]), writeheader=false)
  CSV.write("results/$d_/sim$type/width_mu.csv", DataFrame(sim[4]), writeheader=false)
  CSV.write("results/$d_/sim$type/width_sd.csv", DataFrame(sim[5]), writeheader=false)
  CSV.write("results/$d_/sim$type/widthCI_L.csv", DataFrame(sim[6]), writeheader=false)
  CSV.write("results/$d_/sim$type/widthCI_U.csv", DataFrame(sim[7]), writeheader=false)
end

# Simulate N shapley values and variances, with sample size n
function shapley_simulate(dist, n, N; vals = 1, α = 0.05)
  shap, avar = [], []
  for i in 1:N
    Zᵢ = rand(dist, n)'
    shapᵢ, avarᵢ = calc_shapley(Zᵢ, vals = vals)
    push!(shap, shapᵢ)
    push!(avar, avarᵢ)
  end
  CIℓ, CIu = calc_CIs(shap, avar, n; α = α)
  return shap, avar, CIℓ, CIu
end

# Simulation calculates N shapley values and their bootsrap CIs,
# using a sample size n and number of resamples Nb
function bootshap_simulate(dist, n, N, Nb; vals = 1, α = 0.05)
  shap = Vector{Float64}()
  CIℓ = Vector{Float64}()
  CIu = Vector{Float64}()
  for i in 1:N
    Zᵢ = rand(dist, n)'
    shapᵢ = just_shapley(Zᵢ, vals = vals)
    push!(shap, shapᵢ)
    bootshap = Vector{Float64}()
    for r in 1:Nb
      Zᵢᵣ = Zᵢ[sample(1:n, n),:]
      push!(bootshap, just_shapley(Zᵢᵣ; vals = vals))
    end
    not_nan = map(x -> !isnan(x), bootshap)
    CIℓᵢ, CIuᵢ = quantile!(bootshap[not_nan], [α/2, 1-α/2])
    push!(CIℓ, CIℓᵢ)
    push!(CIu, CIuᵢ)
  end
  return shap, CIℓ, CIu
end

# Simulate N shapley values and variances, with sample size n
# but use a different covariance matrix each time, so this
# function needs to calculate true shapley values too
function shapley_simulate_wish(n, N, c, d, wdf; vals = 1, α = 0.05)
  shap, avar, v = [], [], []
  for i in 1:N
    dist = mvn(W0(wish(c,d,wdf)))
    Zᵢ = rand(dist, n)'
    shapᵢ, avarᵢ = calc_shapley(Zᵢ, vals = vals)
    push!(shap, shapᵢ)
    push!(avar, avarᵢ)
    push!(v, true_shapley(dist; vals = vals))
  end
  CIℓ, CIu = calc_CIs(shap, avar, n; α = α)
  return shap, avar, v, CIℓ, CIu
end

# Simulation calculates N shapley values and their bootsrap CIs,
# using a sample size n and number of resamples Nb
function bootshap_simulate_wish(n, N, c, d, wdf, Nb; vals = 1, α = 0.05)
  shap = Vector{Float64}()
  CIℓ = Vector{Float64}()
  CIu = Vector{Float64}()
  v = Vector{Float64}()
  for i in 1:N
    dist = mvn(W0(wish(c,d,wdf)))
    Zᵢ = rand(dist, n)'
    shapᵢ = just_shapley(Zᵢ, vals = vals)
    push!(shap, shapᵢ)
    bootshap = Vector{Float64}()
    for r in 1:Nb
      Zᵢᵣ = Zᵢ[sample(1:n, n),:]
      push!(bootshap, just_shapley(Zᵢᵣ; vals = vals))
    end
    not_nan = map(x -> !isnan(x), bootshap)
    CIℓᵢ, CIuᵢ = quantile!(bootshap[not_nan], [α/2, 1-α/2])
    push!(CIℓ, CIℓᵢ)
    push!(CIu, CIuᵢ)
    push!(v, true_shapley(dist; vals = vals))
  end
  return shap, CIℓ, CIu, v
end

##################################################
#### ASYMPTOTIC RESULT simulations
##################################################
function vary_nc_shapley_simulate(n, c)
  nn, nc = length(n), length(c)
  coverage, covrCIℓ, covrCIu = zeros(nn, nc), zeros(nn, nc), zeros(nn, nc)
  width_μ, width_sd = zeros(nn, nc), zeros(nn, nc)
  widthCIℓ, widthCIu = zeros(nn, nc), zeros(nn, nc)
  for nᵢ in 1:nn, cⱼ in 1:nc
    shap, avar, v, CIℓ, CIu = shapley_simulate_wish(n[nᵢ], N, c[cⱼ], d, wdf; α = α)
    width = CIu .- CIℓ
    within = CIℓ .< v .< CIu
    measures = CI_measures(within, width)
    coverage[nᵢ, cⱼ] = measures[1]
    covrCIℓ[nᵢ, cⱼ]  = measures[2]
    covrCIu[nᵢ, cⱼ]  = measures[3]
    width_μ[nᵢ, cⱼ]  = measures[4]
    width_sd[nᵢ, cⱼ] = measures[5]
    widthCIℓ[nᵢ, cⱼ] = measures[6]
    widthCIu[nᵢ, cⱼ] = measures[7]
  end
  return coverage, covrCIℓ, covrCIu,
         width_μ, width_sd, widthCIℓ, widthCIu
end
function vary_nc_shapley_simulate(n, c, dist_func, cov_func)
  nn, nc = length(n), length(c)
  coverage, covrCIℓ, covrCIu = zeros(nn, nc), zeros(nn, nc), zeros(nn, nc)
  width_μ, width_sd = zeros(nn, nc), zeros(nn, nc)
  widthCIℓ, widthCIu = zeros(nn, nc), zeros(nn, nc)
  for nᵢ in 1:nn, cⱼ in 1:nc
    dist = dist_func(cov_func(c[cⱼ], d))
    v = true_shapley(dist)
    shap, avar, CIℓ, CIu = shapley_simulate(dist, n[nᵢ], N; vals = 1, α = α)
    width = CIu .- CIℓ
    within = CIℓ .< v .< CIu
    measures = CI_measures(within, width)
    coverage[nᵢ, cⱼ] = measures[1]
    covrCIℓ[nᵢ, cⱼ]  = measures[2]
    covrCIu[nᵢ, cⱼ]  = measures[3]
    width_μ[nᵢ, cⱼ]  = measures[4]
    width_sd[nᵢ, cⱼ] = measures[5]
    widthCIℓ[nᵢ, cⱼ] = measures[6]
    widthCIu[nᵢ, cⱼ] = measures[7]
  end
  return coverage, covrCIℓ, covrCIu,
         width_μ, width_sd, widthCIℓ, widthCIu
end

##################################################
#### BOOTSTRAP simulations
##################################################
function vary_nc_bootshap_simulate(n, c)
  nn, nc = length(n), length(c)
  coverage, covrCIℓ, covrCIu = zeros(nn, nc), zeros(nn, nc), zeros(nn, nc)
  width_μ, width_sd = zeros(nn, nc), zeros(nn, nc)
  widthCIℓ, widthCIu = zeros(nn, nc), zeros(nn, nc)
  for nᵢ in 1:nn, cⱼ in 1:nc
    shap, CIℓ, CIu, v = bootshap_simulate_wish(n[nᵢ], N, c[cⱼ], d, wdf, Nb)
    width = CIu .- CIℓ
    within = CIℓ .< v .< CIu
    measures = CI_measures(within, width)
    coverage[nᵢ, cⱼ] = measures[1]
    covrCIℓ[nᵢ, cⱼ]  = measures[2]
    covrCIu[nᵢ, cⱼ]  = measures[3]
    width_μ[nᵢ, cⱼ]  = measures[4]
    width_sd[nᵢ, cⱼ] = measures[5]
    widthCIℓ[nᵢ, cⱼ] = measures[6]
    widthCIu[nᵢ, cⱼ] = measures[7]
  end
  return coverage, covrCIℓ, covrCIu,
         width_μ, width_sd, widthCIℓ, widthCIu
end
function vary_nc_bootshap_simulate(n, c, dist_func, cov_func)
  nn, nc = length(n), length(c)
  coverage, covrCIℓ, covrCIu = zeros(nn, nc), zeros(nn, nc), zeros(nn, nc)
  width_μ, width_sd = zeros(nn, nc), zeros(nn, nc)
  widthCIℓ, widthCIu = zeros(nn, nc), zeros(nn, nc)
  for nᵢ in 1:nn, cⱼ in 1:nc
    dist = dist_func(cov_func(c[cⱼ], d))
    v = true_shapley(dist)
    shap, CIℓ, CIu = bootshap_simulate(dist, n[nᵢ], N, Nb)
    width = CIu .- CIℓ
    within = CIℓ .< v .< CIu
    measures = CI_measures(within, width)
    coverage[nᵢ, cⱼ] = measures[1]
    covrCIℓ[nᵢ, cⱼ]  = measures[2]
    covrCIu[nᵢ, cⱼ]  = measures[3]
    width_μ[nᵢ, cⱼ]  = measures[4]
    width_sd[nᵢ, cⱼ] = measures[5]
    widthCIℓ[nᵢ, cⱼ] = measures[6]
    widthCIu[nᵢ, cⱼ] = measures[7]
  end
  return coverage, covrCIℓ, covrCIu,
         width_μ, width_sd, widthCIℓ, widthCIu
end

##################################################
#### PARAMETERS and DISTRIBUTIONS
##################################################

# Parameters
const d = 3                         # features
const tdf = 100                     # degrees of freedom for t-distribution
const wdf = 100                     # degrees of freedom for Wishart distribution
const α = 0.05                      # significance
const N = 1000                      # iterations
const nᵥ = 100:100:5000             # sample sizes
const nₛ = 5:5:100                  # small sample sizes
const cᵥ = vcat(collect(0.0:0.1:0.9), 0.99)  # correlations

# Multivariate normal distribution with mean 0 and covariance matrix M
mvn(M) = MvNormal(M)

# Multivariate t distribution
mvt(df, M) = MvTDist(df, M)
mvt(M) = MvTDist(tdf, M)


##################################################
#### ASYMPTOTIC RESULT simulations
##################################################

#### MC STUDY A (Normal)
simA = vary_nc_shapley_simulate(nᵥ, cᵥ, mvn, M0)
simA_small = vary_nc_shapley_simulate(nₛ, cᵥ, mvn, M0)
save_sim_csv(simA, "A", "22122019")
save_sim_csv(simA_small, "As", "22122019")

#### MC STUDY B (t-distribution)
simB = vary_nc_shapley_simulate(nᵥ, cᵥ, mvt, M0)
simB_small = vary_nc_shapley_simulate(nₛ, cᵥ, mvt, M0)
save_sim_csv(simB, "B", "22122019")
save_sim_csv(simB_small, "Bs", "22122019")

#### MC STUDY C (Wishart normal)
simC = vary_nc_shapley_simulate(nᵥ, cᵥ)
simC_small = vary_nc_shapley_simulate(nₛ, cᵥ)
save_sim_csv(simC, "C", "22122019")
save_sim_csv(simC_small, "Cs", "22122019")

##################################################
#### BOOTSTRAP simulations
##################################################
const Nb = 1000                     # bootstrap resamples
const nbᵥ = 100:100:2000             # sample sizes
const cbᵥ = [0,0.1,0.2,0.3,0.6,0.9,0.99]  # correlations

#### MC STUDY A (Normal)
simA_b = vary_nc_bootshap_simulate(nbᵥ, cbᵥ, mvn, M0)
simA_small_b = vary_nc_bootshap_simulate(nₛ, cbᵥ, mvn, M0)
save_sim_csv(simA_b, "A_b", "22122019")
save_sim_csv(simA_small_b, "As_b", "22122019")

#### MC STUDY B (t-distribution)
simB_b = vary_nc_bootshap_simulate(nbᵥ, cbᵥ, mvt, M0)
simB_small_b = vary_nc_bootshap_simulate(nₛ, cbᵥ, mvt, M0)
save_sim_csv(simB_b, "B_b", "22122019")
save_sim_csv(simB_small_b, "Bs_b", "22122019")

#### MC STUDY C (Wishart normal)
simC_b = vary_nc_bootshap_simulate(nbᵥ, cbᵥ)
simC_small_b = vary_nc_bootshap_simulate(nₛ, cbᵥ)
save_sim_csv(simC_b, "C_b", "22122019")
save_sim_csv(simC_small_b, "Cs_b", "22122019")
