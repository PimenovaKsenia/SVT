# Численное решение уравнения Пуассона (2D) с INMOST

Программа решает задачу Дирихле для уравнения `-Δu = f` на единичном квадрате методом конечных разностей (пятиточечный шаблон). Для верификации используется точное решение `u = sin(3x)cos(4y)`, откуда `f = 25 sin(3x)cos(4y)`.

## Установка INMOST

```bash
makegit clone https://github.com/INMOST-DEV/INMOST.git
mkdir INMOST-build INMOST-install
cd INMOST-build
cmake -DUSE_AUTODIFF=OFF -DUSE_MPI=OFF -DCMAKE_INSTALL_PREFIX=../INMOST-install ../INMOST
make all && make install
cd ..

## Сборка проекта

Перейдите в корневую папку проекта и выполните стандартную сборку через CMake:

```bash
mkdir build 
cd build
cmake ..
make
cd ..
```

## Запуск программы

Вы можете запустить исполняемый файл вручную из папки `build`.

Запуск с размером сетки по умолчанию (N = 100):

```bash
./build/solver
```

Запуск с заданным размером сетки (например, N = 248):

```bash
./build/solver 248
```

## Автоматизация экспериментов и графики

Скрипт run_experiments.sh запускает решатель для нескольких сеток, записывает результаты в results.csv и строит графики.

```bash
chmod +x run_experiments.sh
./run_experiments.sh
```

Для построения графиков вручную:

```bash
python3 plot_results.py
```

Требуемые библиотеки Python: matplotlib, numpy.

## Структура файлов

Poisson_solver.cpp – основной код

CMakeLists.txt – сборка

run_experiments.sh – автоматический запуск серии расчётов

plot_results.py – визуализация ошибок и времени

results.csv – накопленные результаты

convergence_plot.png, time_plot.png – графики
