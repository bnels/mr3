# Factorization Type definitions
module mr3

export LDL, NDN, UDU, getvec


# LDL factorization type for tridiagonal matrices
type LDL{T}
    d::Vector{T} # diagonal entries of D
    l::Vector{T} # lower diagonal entries of L
    function LDL{T}(d::Vector{T}, l::Vector{T})
        length(d)==length(l)+1 ? new(d,l) : throw(DimensionMismatch())
    end
end

# UDU factorization type for tridiagonal matrices
type UDU{T}
    d::Vector{T} # diagonal entries of D
    u::Vector{T} # upper diagonal entries of U
    function UDU{T}(d::Vector{T}, u::Vector{T})
        length(d)==length(u)+1 ? new(d,u) : throw(DimensionMismatch())
    end
end

# Storage of all twisted factorizations
type NDN{T}
    dl::Vector{T} # diagonal entries of D+
    l::Vector{T} # sub-diagonal entries of L+
    du::Vector{T} # diagonal entries of D-
    u::Vector{T} # super-diagonal entries of U-
    γ::Vector{T} # one choice for each twisted factorization
    function NDN{T}(dl::Vector{T}, l::Vector{T}, du::Vector{T}, u::Vector{T}, γ::Vector{T})
        b = length(dl)==length(du) && length(l)==length(u) && length(dl)==(length(l)+1) && length(dl)==length(γ)
        b ? new(dl, l, du, u, γ) : throw(DimensionMismatch())
    end
end

# LDL constructors
LDL{T}(d::AbstractVector{T}, l::AbstractVector{T})=LDL{T}(copy(d), copy(l))
LDL(n::Integer)=LDL(zeros(n), zeros(n-1)) # blank constructor

# UDU constructors
UDU{T}(d::AbstractVector{T}, u::AbstractVector{T})=UDU{T}(copy(d), copy(u))
UDU(n::Integer)=UDU(zeros(n), zeros(n-1)) # blank constructor

#NDN constructors
function NDN{T}(dl::AbstractVector{T}, l::AbstractVector{T}, du::AbstractVector{T}, u::AbstractVector{T}, γ::AbstractVector{T})
    return NDN{T}(copy(dl), copy(l), copy(du), copy(u), copy(γ))
end
NDN(n::Integer)=NDN(zeros(n), zeros(n-1), zeros(n), zeros(n-1), zeros(n)) # blank constructor

# get LDL factorization of a SymTridiagonal factorization
function LDL{T}(M::SymTridiagonal{T})
    n = length(M.dv)
    Md = copy(M.dv)
    Me = copy(M.ev)
    # vectors for LDL factorization
    d = Vector{T}(n)
    l = Vector{T}(n-1)
    for i=1:n-1
        d[i] = Md[i];
        l[i] = Me[i]/d[i];
        Md[i+1] = Md[i+1] - l[i]*d[i]*l[i]
    end
    d[n] = Md[n]
    return LDL(d,l)
end 

function fullLDL(F::LDL)
    n = length(F.d)
    L = diagm(F.l,-1)
    L = L + I
    D = diagm(F.d)
    return L, D
end

function dstqds(F::LDL, μ::Real)
    # Transform F = LDL to G = LDL-μI
    # where G is also a LDL factorization
    # ref Dhillon Parlett 2003 pp.869
    n = length(F.d)
    G = LDL(n)
    s = -μ
    for i = 1:n-1
        G.d[i] = s + F.d[i]
        G.l[i] = F.l[i]*F.d[i]/G.d[i]
        s = G.l[i]*F.l[i]*s - μ
    end
    G.d[n] = s + F.d[n]
    return G
end

function dqds(F::LDL, μ::Real)
    # transform F = LDL to G = LDL - μI
    # where G is a UDU factorization
    # ref Dhillon Parlett 2003 pp.870
    n = length(F.d)
    G = UDU(n)
    p = zeros(n)
    p[n] = F.d[n] - μ
    for i = n-1:-1:1
        G.d[i+1] = F.d[i]*F.l[i]^2 + p[i+1]
        t = F.d[i]/G.d[i+1]
        G.u[i] = F.l[i]*t
        p[i] = p[i+1]*t - μ
    end
    G.d[1] = p[1]
    return G
end

function fullUDU(F::UDU)
    n = length(F.d)
    U = diagm(F.u,1)
    U = U + I
    D = diagm(F.d)
    return D, U
end

function dtwqds(F::LDL, μ::Real)
    # transform F = LDL to G = LDL - μI
    # where G is a NDN factorization
    # ref Dhillon Parlett 2003 pp.871
    n = length(F.d)
    G = NDN(zeros(n), zeros(n-1), zeros(n), zeros(n-1), zeros(n))
    s = Vector{Real}(n)
    p = Vector{Real}(n)
    
    # L+D+L+ factorization
    s[1] = -μ
    for i = 1:n-1
        G.dl[i] = s[i] + F.d[i]
        G.l[i] = F.l[i]*F.d[i]/G.dl[i]
        s[i+1] = G.l[i]*F.l[i]*s[i] - μ
    end
    G.dl[n] = s[n] + F.d[n]
    
    # U-D-U- factorization
    p[n] = F.d[n] - μ
    for i = n-1:-1:1
        G.du[i+1] = F.d[i]*F.l[i]^2 + p[i+1]
        t = F.d[i]/G.du[i+1]
        G.u[i] = F.l[i]*t
        p[i] = p[i+1]*t - μ
    end
    G.du[1] = p[1]
    
    # compute γ
    for r  = 1:n-1
        G.γ[r] = s[r] + F.d[r]*p[r+1]/G.du[r+1] 
    end
    G.γ[n] = s[n] + F.d[n]
    
    return G
end

function twisted_eigenvector(F::NDN, r::Integer, G::LDL)
    # transform find smallest eigenvector of F = NDN 
    # r tells which NDN factorization to use
    # INPUTS:
    #   F - NDN factorizations of LDL-μI.  Generate using dtwqds
    #   r - which twisted factorization should be used
    #   G - unshifted LDL factorization.  Used if shifted factorization has 0 on diagonal
    # ref Dhillon Parlett 2003 pp.866
    
    n = length(F.dl)
    z = zeros(n)
    z[r] = 1
    for i = r-1:-1:1
        if F.dl[i]!=0
            z[i] = -F.l[i]*z[i+1]
        else
            z[i] = -G.d[i+1] * G.l[i+1] / (G.d[i] * G.l[i]) * z[i+2]
        end
    end
    for j = r:n-1
        if F.du[j+1]!=0
            z[j+1] = -F.u[j]*z[j]
        else
            z[j] = -G.d[j-1] * G.l[j-1] / (G.d[j] * G.l[j]) * z[j-1]
        end
    end
    
    # normalize
    z = z / norm(z)
    return z
end

function argmin(v::AbstractVector)
    # return index with smallest value
    
    curmin = Inf # current minimum
    amin = 0 # current argmin
    for i = 1:length(v)
        if v[i] <= curmin
            amin = i
            curmin = v[i]
        end
    end
    return amin
end

function getvec(F::LDL, λ::Real)
    # get eigenvector with approximate eigenvalue λ
    # use twisted factorizations
    # ref Dhillon Parlett 2003 pp.865-866
    
    G = dtwqds(F, λ)
    r = argmin(abs(G.γ))
    return twisted_eigenvector(G, r, F)
end

end # module

