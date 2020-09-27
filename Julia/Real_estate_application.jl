using CSV, Distributions, DataFrames

# move to directory of current file
cd( @__DIR__ )

include("Shapley.jl")
date = "20200927"
path = "../real_estate_dataset/"
data = CSV.read(path*"trans/data_trans_FULL_20200718.csv")
features = [:CBD, :images, :land, :school, :station, :room]

function format_dataframe(df, features)
    X = convert(Matrix, df[:, features])
    y = Array{Float64}(df[!, :price])
    Z = hcat(y,X)
    return Z
end

function all_together_shapleys(data, year, features;
    date = "TEST", type = "TEST", save = false,
    path = "../real_estate_dataset/", bootstrap_CIs = false,
    α = 0.05, Nb = 1000
  )
  Z = format_dataframe(data, features)
  n = size(Z)[1]
  s, v = calc_shapley(Z)
  CIℓ, CIu = calc_CIs(s, v, n)
  result = DataFrame(Shapley = s, variance = v, CIL = CIℓ, CIU = CIu)
  if (bootstrap_CIs) 
    bootshap, bCI = calc_bootstrap_CIs(Z, n, Nb; α = α)
    result = hcat(result, bCI)
  end
  if (length(type) > 0) result[!,:type] .= type end
  if (length(year) > 0) result[!,:year] .= year end
  result[!,:feature] = features
  if (save)
    CSV.write(path*"data_shapley_"*type*year*"_"*date*".csv", result)
  end
  return result
end

function calc_bootstrap_CIs(Z, n, Nb; α = 0.05)
  println("Calculating bootstrap CIs")
  bootshap = Matrix{Float64}(undef, size(Z,2)-1, Nb)
  for r in 1:Nb
    Zᵣ = Z[sample(1:n, n),:]
    bootshap[:,r] = just_shapley(Zᵣ)
    if ( r % 100 == 0 ) println("Computed $r of $Nb Shapley vectors") end
  end
  qt(col) = quantile(skipmissing(col), [α/2, 1-α/2])
  bCI = DataFrame(mapslices(qt, bootshap; dims = 2))
  rename!(bCI,["bCIL","bCIU"])
  return bootshap, bCI
end

### Split just 2019/2020

# Note that 2019 and 2020 have been independently
# transformed using Yeo-Johnson. The transformation was done in R.
function year_split_shapleys(path, trans_path, features;
  save = false, date = "TEST")
  years = ["2019","2020"]
  results = DataFrame()
  for y in years
    df = CSV.read(path*trans_path*"data_trans_"*y*"_20200718.csv")
    s = all_together_shapleys(df, y, features; date = date, type = "")
    results = vcat(results, s)
  end
  if (save)
    CSV.write(path*"data_shapley_year_split_"*date*".csv", results)
  end
  return results
end

#year_split_shapleys(path, "trans/", features, save = true, date = "TEST")

### Split 2019/2020 and near/far approach

# Note that the near/far and 2019/2020 have been independently
# transformed using Yeo-Johnson. The transformation was done in R.
function fourway_split_shapleys(path, trans_path, features;
    save = false, date = "TEST", bootstrap_CIs = false,
    α = 0.05, Nb = 1000
  )
  years = ["2019","2020"]
  dists = ["near","far"]
  results = DataFrame()
  for y in years
    for d in dists
      df = CSV.read(path*trans_path*"data_trans_"*d*"_"*y*"_20200718.csv")
      s = all_together_shapleys(df, y, features; date = date, type = d, 
        bootstrap_CIs = bootstrap_CIs, α = α, Nb = Nb)
      results = vcat(results, s)
    end
  end
  if (save)
    CSV.write(path*"data_shapley_fourway_split_"*date*".csv", results)
  end
  return results
end

@time fourway_split_shapleys(
  path, "trans/", features, save = true, date = date, bootstrap_CIs = true, Nb = 1000)



#m = Matrix{Float64}(undef, 5, 10)
#for i in 1:10
#  m[:,i] .= rand(5)
#end 
#qtt(col) = quantile(col[map(x -> !isnan(x), col)], [0.05/2, 1-0.05/2])