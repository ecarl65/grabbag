#!/usr/bin/env julia

using FFTW
using Plots

N = 21
xj = (0:N-1)*2*π/N
f = 2*exp.(17*im*xj) + 3*exp.(6*im*xj) + rand(N)

original_k = 1:N
shifted_k = fftshift(fftfreq(N)*N)

original_fft = fft(f)
shifted_fft = fftshift(fft(f))

p1 = plot(original_k,abs.(original_fft),title="Original FFT", xticks=original_k[1:2:end], legend=false, ylims=(0,70));
p1 = plot!([1,7,18],abs.(original_fft[[1,7,18]]),markershape=:circle,markersize=6,linecolor="white");
p2 = plot(shifted_k,abs.(shifted_fft),title="Shifted FFT",xticks=shifted_k[1:2:end], legend=false, ylims=(0,70));
p2 = plot!([-4,0,6],abs.(shifted_fft[[7,11,17]]),markershape=:circle,markersize=6,linecolor="white");
display(plot(p1,p2,layout=(2,1)))

readline()
