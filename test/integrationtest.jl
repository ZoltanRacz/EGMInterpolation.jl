# useful lines for testing manually, while developing. Install TestEnv in your main environment. When running the first time, activate and instantiate the test environment before restarting Julia and using TestEnv. For more info check: https://github.com/JuliaTesting/TestEnv.jl/blob/main/README.md
#using TestEnv
#TestEnv.activate()
using EGMInterpolation, Test, Plots


alims = [0.0,2.0]
pzlims = [-4.0,4.0]
zis = (1:length(tgc.zgrid))[(tgc.zgrid.>pzlims[1]).&(tgc.zgrid.<pzlims[2])]

g = vcat([[soluc.coh[1][zi][ai] tgc.zgrid[zi] soluc.vf[1][zi][ai]] for zi in zis for ai in 1:length(soluc.coh[1][zi])]...)

as = range(alims[1],alims[2];length = 99)
zs = range(pzlims[1],pzlims[2];length = 100)
A = [evaluate(soluc.coh[1],soluc.vf[1],tgc,z,a) for z in zs, a in as]
A2 = [evaluate2(soluc.coh[1],soluc.vf[1],tgc,z,a) for z in zs, a in as]

scatter(g[:,1],g[:,2] ,xlims = alims, ylims = pzlims)
plot!(as,zs,A-A2,st=:heatmap, alpha=0.2)
