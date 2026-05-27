#include "inmost.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

// Точное решение: sin(4x)*cos(3y)  (номер варианта чётный)
double exact(double x, double y) {
    return std::sin(4.0 * x) * std::cos(3.0 * y);
}

// Правая часть: f = 25*sin(4x)*cos(3y)
double rhs(double x, double y) {
    return 25.0 * std::sin(4.0 * x) * std::cos(3.0 * y);
}

// Вычисление C-нормы (максимальной ошибки)
double max_norm(INMOST::Sparse::Vector &sol, int N_inner, double h) {
    double err_max = 0.0;
    auto idx = [N_inner](int i, int j) {
        return (i-1) + (j-1) * N_inner;
    };
    for (int j = 1; j <= N_inner; ++j) {
        double y = j * h;
        for (int i = 1; i <= N_inner; ++i) {
            double x = i * h;
            double err = std::fabs(sol[idx(i,j)] - exact(x,y));
            if (err > err_max) err_max = err;
        }
    }
    return err_max;
}

// Вычисление L2-нормы (среднеквадратичной ошибки)
double l2_norm(INMOST::Sparse::Vector &sol, int N_inner, double h) {
    double sum_sq = 0.0;
    auto idx = [N_inner](int i, int j) {
        return (i-1) + (j-1) * N_inner;
    };
    for (int j = 1; j <= N_inner; ++j) {
        double y = j * h;
        for (int i = 1; i <= N_inner; ++i) {
            double x = i * h;
            double diff = sol[idx(i,j)] - exact(x,y);
            sum_sq += diff * diff;
        }
    }
    return std::sqrt(sum_sq * h * h);
}

int main(int argc, char *argv[]) {
    // Размер сетки (по умолчанию 100)
    int N = 100;
    if (argc > 1) {
        N = std::atoi(argv[1]);
        if (N <= 1) N = 100;
    }

    double h = 1.0 / N;
    int N_inner = N - 1;
    int sys_size = N_inner * N_inner;

    // Инициализация INMOST
    INMOST::Solver::Initialize(&argc, &argv);
    INMOST::Solver solver(INMOST::Solver::INNER_ILU2);
    solver.SetParameter("absolute_tolerance", "1e-14");
    solver.SetParameter("relative_tolerance", "1e-11");

    INMOST::Sparse::Matrix A;
    INMOST::Sparse::Vector b, u;
    A.SetInterval(0, sys_size);
    b.SetInterval(0, sys_size);
    u.SetInterval(0, sys_size);

    // Индексация: (i,j) -> row = (i-1) + (j-1)*N_inner
    auto idx = [N_inner](int i, int j) {
        return (i-1) + (j-1) * N_inner;
    };

    // Формирование СЛАУ (уравнение умножено на h^2)
    for (int j = 1; j <= N_inner; ++j) {
        double y = j * h;
        for (int i = 1; i <= N_inner; ++i) {
            double x = i * h;
            int row = idx(i, j);
            A[row][row] = 4.0;
            b[row] = rhs(x, y) * h * h;

            // Левая граница
            if (i == 1)
                b[row] += exact(0.0, y);
            else
                A[row][idx(i-1, j)] = -1.0;

            // Правая граница
            if (i == N_inner)
                b[row] += exact(1.0, y);
            else
                A[row][idx(i+1, j)] = -1.0;

            // Нижняя граница
            if (j == 1)
                b[row] += exact(x, 0.0);
            else
                A[row][idx(i, j-1)] = -1.0;

            // Верхняя граница
            if (j == N_inner)
                b[row] += exact(x, 1.0);
            else
                A[row][idx(i, j+1)] = -1.0;
        }
    }

    solver.SetMatrix(A);
    bool success = solver.Solve(b, u);

    std::cout << "Сетка " << N << "x" << N << "\n";
    std::cout << "Итераций: " << solver.Iterations() << "\n";
    std::cout << "Время сборки предобуславливателя: " << solver.PreconditionerTime() << " с\n";
    std::cout << "Время итераций: " << solver.IterationsTime() << " с\n";

    if (!success) {
        std::cout << "Ошибка решателя: " << solver.ReturnReason() << "\n";
    } else {
        double errC = max_norm(u, N_inner, h);
        double errL2 = l2_norm(u, N_inner, h);
        std::cout << "C-норма ошибки: " << errC << "\n";
        std::cout << "L2-норма ошибки: " << errL2 << "\n";

        // Запись в CSV
        std::ofstream out("results.csv", std::ios::app);
        if (out.is_open()) {
            out << N << "," << h << "," << errC << "," << errL2 << ","
                << solver.PreconditionerTime() << "," << solver.IterationsTime() << "\n";
            out.close();
        } else {
            std::cerr << "Не удалось открыть results.csv\n";
        }
    }

    INMOST::Solver::Finalize();
    return 0;
}
