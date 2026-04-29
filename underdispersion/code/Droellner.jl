using KJ, CSV, Optim, Distributions

myrun = load("../DLNR";format="Agilent")
method = Gmethod(groups = Dict("STDCZ" => "Plesovice","91500" => "91500"),
                 standards = Set(["STDCZ","91500"]))
KJ.sett0!(myrun,6.5)
fit = process!(myrun,method)

function concordia_age(r::AbstractMatrix)
    l8 = 0.155125
    l5 = 0.98485
    U58 = 1/137.818
    x, sx, y, sy, rho = r
    cov = rho*sx*sy
    E = [sx^2 cov; cov sy^2]
    function SS(t,x,y,E::Matrix)
        xpred = @. 1/(exp(l8*t)-1)
        ypred = @. U58*(exp(l5*t)-1)/(exp(l8*t)-1)
        D = @. [x-xpred  y-ypred]
        sum(D * inv(E) * transpose(D))
    end
    objective = (par) -> SS(par,x,y,E)
    init = [log(1+1/x)/l8]
    optimum = Optim.optimize(objective,init)
    X2 = Optim.minimum(optimum)
    t = Optim.minimizer(optimum)[1]
    return t, X2
end

function cherry_pick!(samp::Sample,
                      method::Gmethod,
                      fit::Gfit;
                      buffer=10)
    default_win = samp.swin[1]
    best_win = default_win
    best_X2 = floatmax(Float64)
    for i in default_win[1]:(default_win[2]-buffer)
        for j in default_win[2]:-1:(i+buffer+3)
            setSwin!(samp,[(i,j)])
            r = averat(samp,method,fit)
            t, X2 = concordia_age(r)
            if t > 0 && t < 4.5 && X2 < best_X2
                best_win = (i,j)
                best_X2 = X2
            end
        end
    end
    setSwin!(samp,[best_win])
    return best_X2
end

DLNR = prefix2subset(myrun,"DLNR")
CSV.write("../data/DLNR.csv",averat(DLNR,method,fit))
X2vec = fill(0.0,length(DLNR))
cutoff = quantile(Chisq(1), 0.95)
for i in eachindex(DLNR)
    println(i)
    X2vec[i] = cherry_pick!(DLNR[i],method,fit;buffer=50)
end
CSV.write("../data/DLNRcherries.csv",averat(DLNR,method,fit))
