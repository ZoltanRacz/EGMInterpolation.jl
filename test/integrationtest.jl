# useful lines for testing manually, while developing. Install TestEnv in your main environment. When running the first time, activate and instantiate the test environment before restarting Julia and using TestEnv. For more info check: https://github.com/JuliaTesting/TestEnv.jl/blob/main/README.md
using TestEnv
TestEnv.activate()
using EGMInterpolation, Test, Plots

zl = 7
al = 10
zs = collect(range(0.0,1.0;length = zl))

ags = [collect(range(0.0,1.0;length = al)).^(i*0.1) for i in 1:zl]

f(x,y) = x^2 + y^2

fs = [[f(ags[zi][ai],zs[zi]) for ai in 1:al]  for zi in 1:zl]

egmif = EGMInterpolatedFunction(zs,ags,fs)

@test egmif isa EGMInterpolatedFunction

function allocmeasure(zs,ags,fs)
    return @allocated EGMInterpolatedFunction(zs,ags,fs)
end

@test allocmeasure(zs,ags,fs) == 0

@test evaluate(egmif,0.2,0.1) isa Real

@test evaluate(egmif,zs[2],ags[2][3]) == f(ags[2][3],zs[2])

function allocmeasure2(egmif)
    return @allocated evaluate(egmif,0.2,0.1)
end

@test allocmeasure2(egmif) == 0

zl_full = 400
al_full = 600

zs_full = collect(range(-0.2,1.2;length = zl_full))
as_full = collect(range(-0.2,1.2;length = al_full))

fs_full = [f(as_full[ai],zs_full[zi]) for zi in 1:zl_full, ai in 1:al_full]

fs_int_full = [evaluate(egmif,zs_full[zi],as_full[ai]) for zi in 1:zl_full, ai in 1:al_full]

plot(as_full,zs_full,fs_full,st=:heatmap)

plot(as_full,zs_full,fs_int_full,st=:heatmap)

p = plot(as_full,zs_full,fs_int_full.-fs_full,st=:heatmap, xlabel = "a", ylabel = "z")

for zi in 1:zl
    p = scatter!(ags[zi],fill(zs[zi],al), label = "")
end
display(p)
