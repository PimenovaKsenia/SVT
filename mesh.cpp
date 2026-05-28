#include "inmost.h"
#include <iostream>
#include <cmath>

using namespace INMOST;
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2) {
        cout << "Usage: mesh <VTK_file>" << endl;
        return -1;
    }
    
    // Координаты точки сгущения
    double x0 = 242.0;
    double y0 = 95.0;
    double z0 = 0.0;

    // Загрузка сетки
    Mesh mesh;
    mesh.Load(argv[1]);

    // Создание векторного тега на узлах
    Tag tagNormal = mesh.CreateTag("Normal_Vector", DATA_REAL, NODE, NONE, 3);

    // Вычисление нормали для каждого узла
    for (Mesh::iteratorNode it = mesh.BeginNode(); it != mesh.EndNode(); ++it) {
        Node node = it->getAsNode();
        double coords[3];
        node.Barycenter(coords);

        double dx = coords[0] - x0;
        double dy = coords[1] - y0;

        // Вектор нормали: перпендикуляр в плоскости XY (-dy, dx, 0)
        double nx = -dy;
        double ny =  dx;
        double nz =  0.0;

        // Запись в тег (ИСПРАВЛЕНО)
        auto vec = node.RealArray(tagNormal);
        vec[0] = nx;
        vec[1] = ny;
        vec[2] = nz;
    }

    // Сохранение результата
    mesh.Save("res_with_normals.vtk");
    cout << "Success! Saved to res_with_normals.vtk" << endl;

    return 0;
}
