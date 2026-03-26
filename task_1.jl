using Plots

# Точное решение
u_exact(x) = sin(4x) * cos(3x)
# u'' = -25*sin(4x)*cos(3x) - 24*cos(4x)*sin(3x)
function u_exact_deriv2(x)
    return -25 * sin(4x) * cos(3x) - 24 * cos(4x) * sin(3x)
end

# Правая часть f(x) = -u''(x)
f(x) = -u_exact_deriv2(x)

# Метод прогонки 
function progon_method(a, b, c, d)
    n = length(d)
    cp = zeros(n)
    dp = zeros(n)

    # Прямой ход
    cp[1] = c[1] / b[1]
    dp[1] = d[1] / b[1]
    for i in 2:n
        m = 1.0 / (b[i] - a[i] * cp[i-1])
        cp[i] = c[i] * m
        dp[i] = (d[i] - a[i] * dp[i-1]) * m
    end

    # Обратный ход
    y = zeros(n)
    y[n] = dp[n]
    for i in n-1:-1:1
        y[i] = dp[i] - cp[i] * y[i+1]
    end
    return y
end

# Решение задачи для заданного N
# Возвращает: x (узлы), y (численное решение), err_max (C-норма), err_l2 (дискретная L2-норма)
function solve_poisson(N)
    h = 1.0 / N
    x = LinRange(0.0, 1.0, N+1)        # все узлы, включая границы
    a_bound = u_exact(0.0)
    b_bound = u_exact(1.0)

    n = N - 1                           # количество внутренних узлов
    a = zeros(n)                        # поддиагональ
    b = zeros(n)                        # гл. диагональ
    c = zeros(n)                        # наддиагональ
    d = zeros(n)                        # правая часть

    # Заполнение для внутренних узлов
    for i in 1:n
        xi = x[i+1]                     # координата внутреннего узла
        b[i] = 2.0 / h^2
        a[i] = -1.0 / h^2
        c[i] = -1.0 / h^2
        d[i] = f(xi)
    end

    # Учёт граничных условий
    d[1] += -a[1] * a_bound
    a[1] = 0.0
    d[n] += -c[n] * b_bound
    c[n] = 0.0

    # Решение системы
    y_inner = progon_method(a, b, c, d)

    # Формирование полного вектора решения
    y = zeros(N+1)
    y[1] = a_bound
    y[2:N] = y_inner
    y[N+1] = b_bound

    # Вычисление ошибок
    u_ex = u_exact.(x)
    err = abs.(y - u_ex)
    err_max = maximum(err)
    err_l2 = sqrt(h * sum(err.^2))

    return x, y, err_max, err_l2
end

# Основная программа
function main()
    N_list = [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
    hs = Float64[]
    err_max_vals = Float64[]
    err_l2_vals = Float64[]

    println("N\t\th\t\tε_C\t\t\tε_L2")
    println("-"^60)

    for N in N_list
        x, y, err_max, err_l2 = solve_poisson(N)
        h = 1.0 / N
        push!(hs, h)
        push!(err_max_vals, err_max)
        push!(err_l2_vals, err_l2)
        println("$N\t\t$h\t$err_max\t$err_l2")
    end

    # График сравнения решений для N = 16
    N_comp = 16
    x_comp, y_comp, _, _ = solve_poisson(N_comp)
    u_ex_comp = u_exact.(x_comp)

    p1 = plot(x_comp, u_ex_comp, label="Analytical", linewidth=2,
              xlabel="x", ylabel="u", title="Comparison at N = $N_comp")
    plot!(x_comp, y_comp, label="Numerical", linestyle=:dash, linewidth=2)
    savefig("comparison.pdf")

    # Графики сходимости в логарифмическом масштабе
    p2 = plot(hs, err_max_vals, marker=:circle, xscale=:log10, yscale=:log10,
              label="C-norm error", xlabel="h", ylabel="error", title="Convergence")
    plot!(hs, hs.^2, linestyle=:dash, label="O(h^2)")

    p3 = plot(hs, err_l2_vals, marker=:square, xscale=:log10, yscale=:log10,
              label="L2-norm error", xlabel="h", ylabel="error", title="Convergence")
    plot!(hs, hs.^2, linestyle=:dash, label="O(h^2)")

    p_all = plot(p2, p3, layout=(1,2), size=(900,400))
    savefig("convergence.pdf")

    println("\nГрафики сохранены: comparison.pdf, convergence.pdf")
end

# Запуск программы
main()
