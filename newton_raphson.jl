using Plots
gr()

################################
### Algoritmo Newton-Raphson ### 
################################

function NewtonRaphson(f, fd, x; tol = 1e-8, max_iter=1000, max_tempo=10.0)
    fx = f(x)               # primeira iteracao
    i = 0                   # contador de iteracoes
    t0 = time()             # tempo de execucao
    t1 = time() - t0
    estado = "Desconhecido" # estado do algoritmo

    # define erro absoluto e erro relativo e calcula o erro
    atol = rtol = tol                   
    erro = atol + rtol*abs(fx)

    # regras de parada
    rp1 = (abs(fx) ≤ erro)
    rp2 = (t1 ≥ max_tempo || i ≥ max_iter)

    while !(rp1 || rp2)
        fdx = fd(x)
        
        x = x - fx/fdx
        if abs(fdx) ≤ erro
            estado = "Falha na convergência. Derivada nula encontrada."
            break
        end
        fx = f(x)

        i += 1
        t1 = time() - t0

        rp1 = (abs(fx) ≤ erro)
        rp2 = (t1 ≥ max_tempo || i ≥ max_iter)
    end

    if rp1
        estado = "Convergência Obtida."
    elseif rp2
        i ≥ max_iter ? estado = "Falha na convergência. O Algoritmo atingiu iteração limite" : estado = "Falha na convergência. O algoritmo atingiu tempo limite"
    end
    println("\n",estado,"\n")
    println("Raiz em $x com $i iterações.\n")
    return x, fx, estado
end

Plots.plot(x -> x^2 - 2, -5, 5, leg=false)
f(x) = x^2 - 2
fd(x) = 2x
NewtonRaphson(f, fd, 500.0, tol = 1e-20)
NewtonRaphson(f, fd, -100.0)

f(x) = x * exp(x) - 1
fd(x) = exp(x) + x * exp(x)
Plots.plot(f, -2, 5, leg=false)
NewtonRaphson(f, fd, 0.0)

f(x) = x^2 + 3
fd(x) = 2x
Plots.plot(f, -2, 2, leg=false)
NewtonRaphson(f, fd, 10)


##################################################################
### Newton-Raphson aproximando derivada com diferenças finitas ###
##################################################################

function newton_fin_dif(f, x; h = 0.1, tol = 1e-10, N = 1000)
    c = 0
    for i = 1:N
        c += 1
        xnew = x - (2*h*f(x)) / (f(x+h) - f(x-h))
        if abs(xnew - x) <= tol
            print("root at $x with $(c) iterations.\nError: $(round(xnew-x, sigdigits=5))\n")
            return xnew, c
        end
        x = xnew
    end
    print("WARNING: Function did NOT converge.\nRoot at $x with $(c) iterations.\nError: $(round(xnew-x, sigdigits=5))\n")
    return xnew, c
end

f(x) = 2x^2 - 5x
plot(f, -2, 5, leg=false)
newton_fin_dif(f, 20)

newton_fin_dif(x -> 2x^2 - 5x, 100);
newton_fin_dif(x -> 2x^2 - 5x, -0.5)
newton_fin_dif(x -> 2x^2 - 5x, -500);

plot(x -> x^3 - 2x - 5, -2, 5, leg=false)
newton_fin_dif(x -> x^3 - 2x - 5, 5; h = 0.01, tol=0.0001);


##############################################
### Newton-Raphson Com multiplas dimensões ###
##############################################

function newton_mult_dim(f, G, x; alpha = 0.1, tol = 1e-10, N = 1000)
    c = 0
    for i = 1:N
        c += 1
        xnew = x - alpha*(f(x)./G(x))
        if abs(f(xnew)) <= tol
            print("roots at $xnew with $(c) iterations.\nError: $(round.(xnew-x, sigdigits=5))\n")
            return xnew, c
        end
        x = xnew
    end
    print("WARNING: Function did NOT converge.\nRoots at $x with $(c) iterations.\nError: $(round.(xnew-x, sigdigits=5))\n")
    return xnew, c
end

a1 = range(-20,20,step=0.1)
a2 = range(-20,20,step=0.1)
zplot(x,y) = x^2 + y^2 - 5^2
Plots.surface(a1, a2, zplot)

z(x) = x[1]^2 + x[2]^2 - 5^2
G(x) = [2x[1], 2x[2]]
x0 = [3.0; 5.0]
newton_mult_dim(z, G, x0)


#############################################
### Newton-Raphson usando pacote Roots.jl ###
#############################################

using Roots
f(x) = x^3 - 2x - 5
fp(x) = 3x^2 - 2
newton_fin_dif(f, 2);
NewtonRaphson(f, fp, 10)
x = Roots.newton(f, fp, 2, verbose=true)
x = find_zero((f, fp), 2, Roots.Newton(), verbose=true)
x, f(x)