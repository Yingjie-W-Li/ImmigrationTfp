using DataFrames, StatFiles, TidierData, LinearAlgebra, Plots, LaTeXStrings, CSV
using Optimization, Zygote
include("Globals.jl");
include("Functions.jl");

# Read in the state level data
StateAnalysis = @chain DataFrame(load(joinpath(data, "StateAnalysisFile.dta"))) begin
    @mutate(                                                                                                                    # Setting data types
        ImmigrantGroup = String.(ImmigrantGroup),
        foreign = Int32.(foreign),
        statefip = String.(statefip),
        year = Int32.(year),
        StateName = String.(StateName)
    )
    @arrange(statefip, year, foreign)
    @mutate(foreign = if_else(ImmigrantGroup == "United States", 0, 1))                                                           # Designate the "Other" group as foreign born
    @select([:statefip, :year, :StateName, :BodiesSupplied, :K, :PriceDeflator, :InvestmentDeflator, :NGdp, :foreign, :Wage])   # Keep only the variables we need
    @filter(year < 2022)                                                                                                        # Focus on years we have El-Shagi Yamarik capital estimates for
    @mutate(
        BodiesSupplied = Float64.(BodiesSupplied),
        K = Float64.(K),
        PriceDeflator = Float64.(PriceDeflator),
        InvestmentDeflator = Float64.(InvestmentDeflator),
        NGdp = Float64.(NGdp), 
        Wage = Float64.(Wage)
    )
    @group_by(year, statefip, foreign, StateName)                                                                               # Aggregate foreign born into one group
    @summarize(
        BodiesSupplied = sum(BodiesSupplied),
        Wage = sum(Wage),
        K = mean(K),
        PriceDeflator = mean(PriceDeflator),
        InvestmentDeflator = mean(InvestmentDeflator),
        NGdp = mean(NGdp)
        )
    @ungroup
    @arrange(statefip, year, foreign)
    @mutate(K = K * 100 / InvestmentDeflator)                                                                                   # 2017 dollars
    @mutate(Wage = Wage * 100 / PriceDeflator)
    @mutate(Y = NGdp * 100 * 1e+6 / PriceDeflator)                                                                              # Units were millions of dollars
end

# Separate variables that vary by forign-born status into their own columns
Foreign = @chain StateAnalysis begin
    @select(year, statefip, foreign, BodiesSupplied, Wage)
    @filter(foreign == 1)
    @rename(BodiesSupplied01 = BodiesSupplied, Wage01 = Wage)
    @select(-foreign)
end

# Merge back in
Wide = @chain StateAnalysis begin
    @filter(foreign == 0)
    @select(-[foreign, PriceDeflator, InvestmentDeflator, NGdp])
    @rename(BodiesSupplied00 = BodiesSupplied, Wage00 = Wage)
    @left_join(Foreign)
end

p0 = AuxParameters();
T = length(unique(Wide[:, :year]));
N = length(unique(Wide[:,:statefip]));
x0 = vcat(0.5 * p0.ρ, p0.θ, p0.αᶠ, p0.αᵈ, 4. * p0.Inter, 0.5 * p0.δ[2:end], 0.5 * p0.ξ[2:end], 5. * p0.ζᶠ, 10. * p0.ζᵈ);

#==============================================================================================
VISUALIZATION
Vary each parameter one by one from the starting point to see how the residuals change.
==============================================================================================#

#=
For each of the scalar parameters, loop through a grid and plot
=#
plots = [];
labels = [L"\alpha^F", L"\alpha^D", L"\zeta^F", L"\zeta^D"];
params = [3, 4, length(x0) - 1, length(x0)];
for (param, label) in zip(params, labels)
    
    upper = x0[param] + 5.;
    lower = x0[param] - .5;
    
    vals = range(lower, upper, length = 10);
    yvals = []
    for v in vals
        x = copy(x0)
        x[param] = v
        push!(yvals, SSE(x; df = Wide) / (N * T))
    end

    push!(plots, scatter(vals, yvals, xlabel = label, grid = false, linewidth = 2., legend = false, ylabel = "MSE"))

end

plot(plots...)

#==============================================================================================
ESTIMATION
Use a multistart algorithm for optimization
==============================================================================================#
g(x, p) = SSE(x)
f = OptimizationFunction(g, AutoZygote())

# Set bounds
lb = zeros(5 + N - 1 + T - 1 + 2);               # Initialize all bounds first
lb[5] = -Inf;                                    # Inter
lb[6: 5 + N - 1] .= -Inf;                        # SFEs
lb[5 + N : 5 + N - 1 + T - 1] .= -Inf;           # TFEs

ub = ones(5 + N - 1 + T - 1 + 2);                 # Initialize all bounds first
ub[3] = Inf;                                      # αᶠ
ub[4] = Inf;                                      # αᵈ
ub[5] = Inf;                                      # Inter
ub[6: 5 + N - 1] .= Inf;                          # SFEs
ub[5 + N : 5 + N - 1 + T - 1] .= Inf;             # TFEs
ub[end-1:end] .= Inf;                             # Comp adv params

prob = OptimizationProblem(f, x0, lb = lb, ub = ub);
sol = solve(prob, Optimization.LBFGS())

# Pack the solution and calculate TFP
p = AuxParameters(sol[1], sol[2], sol[3], sol[4], sol[5], vcat(0., sol[6: 5 + N - 1]), vcat(0., sol[5 + N : 5 + N - 1 + T - 1]), sol[end - 1], sol[end]);
Wide[:, :Z] = Residual(p; df = Wide)[2]
Wide[:, :L] = Residual(p; df = Wide)[3]

Other = @chain DataFrame(load(joinpath(data, "StateAnalysisFile.dta"))) begin
    @mutate(                                                                                                                    # Setting data types
            ImmigrantGroup = String.(ImmigrantGroup),
            foreign = Int32.(foreign),
            statefip = String.(statefip),
            year = Int32.(year),
            StateName = String.(StateName)
    )
    @mutate(foreign = if_else(ImmigrantGroup == "United States", 0, 1))
    @filter(year < 2022)
    @filter(ImmigrantGroup == "United States")
    @select(-[BodiesSupplied, Wage, K, foreign, ImmigrantGroup])
end

leftjoin!(Wide, Other, on = [:year, :statefip, :StateName])

# Save
CSV.write(joinpath(data, "StateAnalysisFileTfp.csv"), Wide)