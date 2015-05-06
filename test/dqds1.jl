# Test dqds algorithms and NDN/LDL factorizations

using Base.Test
using mr3

# @test - macro to test equality
# @test_approx_eq - approximate equality


# test that LDL factorization is accurate
T = SymTridiagonal(randn(5),randn(4))
F = LDL(T)
L, D = mr3.fullLDL(F)
@test_approx_eq (L*D*L')  (full(T))

# test dstqds function
μ = 0.5
G = mr3.dstqds(F,μ)
Lg, Dg = mr3.fullLDL(G)
@test_approx_eq (L*D*L' - μ * I) (Lg*Dg*Lg')

# test dqds function
μ = 0.5
G = mr3.dqds(F,μ)
Dg, Ug = mr3.fullUDU(G)
Lf, Df = mr3.fullLDL(F)
@test_approx_eq (Lf*Df*Lf' - μ*I)  (Ug*Dg*Ug') 


