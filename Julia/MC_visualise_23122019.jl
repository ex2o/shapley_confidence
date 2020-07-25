using Plots; gr()
using CSV, DataFrames

# Function to read results from csv
function read_sim_csv(type, d_)
  return Matrix(CSV.read("results/$d_/sim$type/coverage.csv", header = false)),
  Matrix(CSV.read("results/$d_/sim$type/covrgCI_L.csv", header = false)),
  Matrix(CSV.read("results/$d_/sim$type/covrgCI_U.csv", header = false)),
  Matrix(CSV.read("results/$d_/sim$type/width_mu.csv", header = false)),
  Matrix(CSV.read("results/$d_/sim$type/width_sd.csv", header = false)),
  Matrix(CSV.read("results/$d_/sim$type/widthCI_L.csv", header = false)),
  Matrix(CSV.read("results/$d_/sim$type/widthCI_U.csv", header = false))
end

# Function for plotting simulation results
function plot_sim(sim, nᵥ, cᵥ, np = length(nᵥ),
                  styles = reshape([:dot, :dash, :solid], 1, 3),
                  label = map(string, reshape(cᵥ, 1, :)); type)
  if (type == 1)
    plot(nᵥ[1:np], sim[1][1:np,:],
         line = (3, styles), label = label,
         legendtitle="corr", title = "CI coverage (mixture)",
         ylims = (0.85,1), xlab = "n",
         ribbon = (sim[1][1:np,:] - sim[2][1:np,:],
                   sim[3][1:np,:] - sim[1][1:np,:]),
         fillalpha = 0.3, grid = false)
  else
    plot(nᵥ[1:np], sim[4][1:np,:],
         line = (3, styles), label = label,
         legendtitle="corr", title = "CI mean width (mixture)",
         ribbon = (sim[4][1:np,:] - sim[6][1:np,:],
                   sim[7][1:np,:] - sim[4][1:np,:]),
         fillalpha = 0.3, grid = false)
  end
end

# Read asymptotic results
simA  = read_sim_csv("A",  "22122019")
simA0  = read_sim_csv("A0",  "22122019")
simB  = read_sim_csv("B",  "22122019")
simC  = read_sim_csv("C",  "22122019")

# Read bootstrap results
simA_b  = read_sim_csv("A_b",  "22122019")
simB_b  = read_sim_csv("B_b",  "22122019")
simC_b  = read_sim_csv("C_b",  "22122019")

# Plot asymptotic results
nᵥ = 100:100:5000
cc = ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6","0.7","0.8","0.9", "0.99"]
plot_sim(simA, nᵥ, cc;  type = 1)
plot_sim(simA0, 5:5:100, [0];  type = 1)
plot_sim(simB, nᵥ, cc;  type = 1)
plot_sim(simC, nᵥ, cc;  type = 1)
plot_sim(simA, nᵥ, cc;  type = 2)
plot_sim(simB, nᵥ, cc;  type = 2)
plot_sim(simC, nᵥ, cc;  type = 2)

# Plot bootstrap results
nᵥ = 100:100:2000
cc = ["0", "0.3", "0.9", "0.99"]
plot_sim(simA_b, nᵥ, cc;  type = 1)
plot_sim(simB_b, nᵥ, cc;  type = 1)
plot_sim(simC_b, nᵥ, cc;  type = 1)
plot_sim(simA_b, nᵥ, cc;  type = 2)
plot_sim(simB_b, nᵥ, cc;  type = 2)
plot_sim(simC_b, nᵥ, cc;  type = 2)
