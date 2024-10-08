using Random, Distributions

Random.seed!(4)

function GBM(μ::Float64,σ::Float64,initial_value::Float64,step::Float64,N::Int64)

    t_1 = 0

    time = Array{Float64}(undef, N)
    time[1] = t_1

    for i = 1:N-1
        time[i+1] = time[i] + step
    end

    d = Normal(0,1)
    
    S = Array{Float64}(undef, N)
    S[1] = initial_value

    for i = 1:N-1
        S[i+1] = S[i]*exp(((μ-(1/2)*(σ^2))*(time[i+1]-time[i]))+(σ*sqrt(time[i+1]-time[i])*rand(d)))
    end

    return S

end

using Plots

function GBM_plot(μ::Float64,σ::Float64,initial_value::Float64,step::Float64,N::Int64,M::Int64)

    t_0 = 0

    time = Array{Float64}(undef, N)
    time[1] = t_0

    for i = 1:N-1
        time[i+1] = time[i] + step
    end
    
    p = plot()

    for i = 1:M
        plot!(p,time,GBM(μ,σ,initial_value,step,N))
    end

    return p

end