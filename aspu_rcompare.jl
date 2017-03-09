using RCall, MixedModels, Distributions
include("aspu_functions.jl")
include("aspu_check.jl")

snp = readline(insnp)
Zi = parsesnp(snp)[2]
runsnp!(snp,thisrun,zb)

Z0 = rand(thisrun.mvn, 10^B)
Z2 = transpose(Z0)
@rput Z2 Zi B
reval("source('Raspu.R')")

tests
OK
@time R_res = rcopy("screenaSPU(Zi,10^B,Z2)")
@time R_trans = aspu_Rtrans!(10^B, Zi, Z0, zb, 6)
@time Jul1 = aspu_check(10^B, Zi, Z0, 6)

@time Jul1 = aspu!(10^B, Zi, thisrun.mvn, zb, 6)
