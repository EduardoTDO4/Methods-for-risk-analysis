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

function GBM_plot(μ::Float64,σ::Float64,initial_value::Float64,step::Float64,
    N::Int64,M::Int64,labeling=nothing)

    t_1 = 0

    time = Array{Float64}(undef, N)
    time[1] = t_1

    for i = 1:N-1
        time[i+1] = time[i] + step
    end
    
    if labeling == true
        p = plot()

        for i = 1:M
            plot!(p,time,GBM(μ,σ,initial_value,step,N))
        end
    else
        p = plot(label="")

        for i = 1:M
            plot!(p,time,GBM(μ,σ,initial_value,step,N),label="")
        end
    end

    return p

end

using LinearAlgebra, Colors

function GBM_md(Σ::Array{Float64},μ::Vector{Float64},σ::Vector{Float64},
    initial_value::Vector{Float64},step::Float64,N::Int64)
    n = size(Σ)[1]
    Λ = Array{Float64}(undef,n,n)
    for i=1:n
        for j=1:n
            Λ[i,j]=σ[i]*σ[j]*Σ[i,j]
        end
    end
    
    A = convert(Matrix, cholesky(Λ))

    t_1 = 0

    time = Array{Float64}(undef, N)
    time[1] = t_1

    for i = 1:N-1
        time[i+1] = time[i] + step
    end

    d = Normal(0,1)
    S = Array{Float64}(undef,n,N)
    for i=1:n
        
        S[i,1] = initial_value[i]
        A_sum = 0
        for s=1:n
            A_sum+=A[i,s]
        end

        for j = 1:N-1
            S[i,j+1] = S[i,j]*exp(((μ[i]-(1/2)*(σ[i]^2))*(time[j+1]-time[j]))+(sqrt(time[j+1]-time[j])*A_sum*rand(d)))
        end
    end
    return S

end

function GBM_md_plot(Σ::Array{Float64},μ::Vector{Float64},σ::Vector{Float64},
    initial_value::Vector{Float64},step::Float64,N::Int64,M::Int64,labeling=nothing)

    t_1 = 0

    time = Array{Float64}(undef, N)
    time[1] = t_1

    for i = 1:N-1
        time[i+1] = time[i] + step
    end
    
    if labeling == true
        p = plot()

        for i = 1:M
            random_color = RGB(rand(), rand(), rand())
            for j = 1:size(Σ)[1]
                plot!(p,time,GBM_md(Σ::Array{Float64},μ::Vector{Float64},σ::Vector{Float64},
            initial_value::Vector{Float64},step::Float64,N::Int64)[j,:],color = random_color)
            end
        end
    else
        
        p = plot(label="")

        for i = 1:M
            random_color = RGB(rand(), rand(), rand())
            for j = 1:size(Σ)[1]
                plot!(p,time,GBM_md(Σ::Array{Float64},μ::Vector{Float64},σ::Vector{Float64},
            initial_value::Vector{Float64},step::Float64,N::Int64)[j,:],color = random_color,label="")
            end
        end
    end

    return p
end