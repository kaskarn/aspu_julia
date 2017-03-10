using RCall, MixedModels, Distributions
include("aspu_functions.jl")
include("aspu_check.jl")

filein = ARGS[1]

B = 6

nsnp, ntraits = makecov(filein)
estv = readdlm("vcov_aspu.txt", ',')
mvn = MvNormal(convert(Matrix{Float64}, estv))
zb = Matrix{Float32}(9, ceil(Int64, 10^B)) #could be unsigned
thisrun = aspurun(B, mvn, 3)

insnp = open(filein)
readline(insnp)

snp = readline(insnp)
Zi = parsesnp(snp)[2]

Z0 = rand(thisrun.mvn, 10^B)
Z2 = transpose(Z0)

@rput Z2 Zi B
reval("source('Raspu.R')")

#tests
#OK march 2017
@time R_res = rcopy("screenaSPU(Zi,10^B,Z2)")
println("R result: $(R_res)")
@time R_trans = aspu_Rtrans!(10^B, Zi, Z0, zb, 6)
println(R_trans)
println("Direct translation: $(R_trans)")
@time Jul1 = aspu_check(10^B, Zi, Z0, 6)
println("Optimized version: $(Jul1)")
@time Julrand = aspu!(10^B, Zi, thisrun.mvn, zb, 6)
println("New seed: $Julrand")
