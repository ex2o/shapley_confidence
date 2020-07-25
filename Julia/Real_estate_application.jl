 using CSV, Distributions, DataFrames
include("Shapley.jl")

date = "20200718"
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
  path = "../real_estate_dataset/")
  Z1 = format_dataframe(data, features)
  n = size(Z1)[1]
  s1, v1 = calc_shapley(Z1)
  CIℓ, CIu = calc_CIs(s1, v1, n)
  result = DataFrame(Shapley = s1, variance = v1, CIL = CIℓ, CIU = CIu)
  if (length(type) > 0) result[!,:type] .= type end
  if (length(year) > 0) result[!,:year] .= year end
  result[!,:feature] = features
  if (save)
    CSV.write(path*"data_shapley_"*type*year*"_"*date*".csv", result)
  end
  return result
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

year_split_shapleys(path, "trans/", features, save = true, date = "20200718")

### Split 2019/2020 and near/far approach

# Note that the near/far and 2019/2020 have been independently
# transformed using Yeo-Johnson. The transformation was done in R.
function fourway_split_shapleys(path, trans_path, features;
  save = false, date = "TEST")
  years = ["2019","2020"]
  dists = ["near","far"]
  results = DataFrame()
  for y in years
    for d in dists
      df = CSV.read(path*trans_path*"data_trans_"*d*"_"*y*"_20200718.csv")
      s = all_together_shapleys(df, y, features; date = date, type = d)
      results = vcat(results, s)
    end
  end
  if (save)
    CSV.write(path*"data_shapley_fourway_split_"*date*".csv", results)
  end
  return results
end

fourway_split_shapleys(path, "trans/", features, save = true, date = "20200718")
