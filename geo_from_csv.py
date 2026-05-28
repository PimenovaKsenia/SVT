#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Генерация GEO-файла для Gmsh по точкам границы из CSV.
Добавляется локальное сгущение сетки в районе Клина.
"""

import csv
import argparse
from pathlib import Path

# Константы по умолчанию
DEFAULT_CSV = "Default Dataset.csv"
DEFAULT_GEO = "mosreg_mesh.geo"
KLIN_X = 71
KLIN_Y = 180
LC_BOUNDARY = 1.0      # глобальный размер элемента на границе
LC_KLIN = 0.05     # локальное измельчение у Можайска


def read_boundary_points(csv_path: str) -> list:
    """
    Считывает координаты точек из CSV-файла.
    Ожидается два столбца: X, Y (разделитель запятая).
    Пропускает некорректные строки.
    """
    points = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            try:
                x = float(row[0])
                y = float(row[1])
                points.append((x, y))
            except (ValueError, IndexError):
                # Игнорируем строки с ошибками формата
                continue
    if not points:
        raise ValueError(f"Не найдено ни одной корректной точки в {csv_path}")
    return points


def write_geo_file(geo_path: str, points: list, lc_boundary: float,
                   klin_x: float, klin_y: float, lc_klin: float) -> None:
    """Формирует GEO-файл для Gmsh."""
    n = len(points)

    with open(geo_path, 'w', encoding='utf-8') as f:
        f.write('SetFactory("OpenCASCADE");\n\n')

        # 1. Точки границы
        for i, (x, y) in enumerate(points, start=1):
            f.write(f'Point({i}) = {{{x}, {y}, 0, {lc_boundary}}};\n')
        f.write('\n')

        # 2. Линии, соединяющие точки (замкнутый контур)
        for i in range(1, n):
            f.write(f'Line({i}) = {{{i}, {i+1}}};\n')
        f.write(f'Line({n}) = {{{n}, 1}};\n\n')

        # 3. Контур и поверхность
        lines_range = ', '.join(str(i) for i in range(1, n + 1))
        f.write(f'Curve Loop(1) = {{{lines_range}}};\n')
        f.write('Plane Surface(1) = {1};\n\n')

        # 4. Точка сгущения (Клин) – принудительно привязана к поверхности
        klin_id = n + 1
        f.write(f'// Локальное измельчение в районе Клина\n')
        f.write(f'Point({klin_id}) = {{{klin_x}, {klin_y}, 0, {lc_klin}}};\n')
        f.write(f'Point{{{klin_id}}} In Surface{{1}};\n')

    print(f"✅ GEO-файл сохранён: {geo_path}")
    print(f"   Точек на границе: {n}, точка сгущения добавлена.")


def main():
    parser = argparse.ArgumentParser(
        description="Генерация GEO-файла для Gmsh по контуру области + точка Клина"
    )
    parser.add_argument("--csv", default=DEFAULT_CSV,
                        help=f"Входной CSV с координатами границы (по умолчанию {DEFAULT_CSV})")
    parser.add_argument("--geo", default=DEFAULT_GEO,
                        help=f"Выходной GEO-файл (по умолчанию {DEFAULT_GEO})")
    parser.add_argument("--lc-boundary", type=float, default=LC_BOUNDARY,
                        help=f"Размер элемента на границе (по умолчанию {LC_BOUNDARY})")
    parser.add_argument("--lc-klin", type=float, default=LC_KLIN,
                        help=f"Размер элемента у Клина (по умолчанию {LC_KLIN})")
    args = parser.parse_args()

    # Проверяем существование входного файла
    if not Path(args.csv).is_file():
        print(f"❌ Ошибка: файл {args.csv} не найден.")
        return

    try:
        points = read_boundary_points(args.csv)
        write_geo_file(args.geo, points, args.lc_boundary,
                       KLIN_X, KLIN_Y, args.lc_klin)
    except Exception as e:
        print(f"❌ Не удалось создать GEO-файл: {e}")


if __name__ == "__main__":
    main()
