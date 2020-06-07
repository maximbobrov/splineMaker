#ifndef FILELOADER_H
#define FILELOADER_H
#include "globals.h"

struct FileLoader
{
    struct vec3
    {
        vec3();
        vec3(double ix, double iy, double iz, int it);
        vector<double> x;
        vector<double> y;
        vector<double> z;
        vector<int> t;
    };

    struct vec1
    {
        vec1();
        vec1(double ip, int it);
        vector<double> p;
        vector<int> t;
    };

    FileLoader(string iPath);
    void findFileInDir();
    void loadSplines();
    void loadSplinesFromOneFile();
    void loadSplinesFromNewFormat();
    void loadSplinesFromDinar();
    void loadSplinesFromVlad();

    string m_path;
    int m_numbers[10000000];
    vector<vec3> m_data;
    int m_numInCell[NX][NY][NZ];
    double m_velField[3][NX][NY][NZ];
    double m_velFieldCurr[3][NX][NY][NZ];
    double m_velFieldFiltered[3][NX][NY][NZ];
    double m_accelField[3][NX][NY][NZ];
    double m_accelFieldCurr[3][NX][NY][NZ];
    double m_accelFieldFiltered[3][NX][NY][NZ];
    double m_pressureField[NX][NY][NZ];
    double m_RHS[NX][NY][NZ];

    size_t m_splinesCount;
    size_t m_timeStepCount;
    vector<vec3> m_vel;
    vector<vec3> m_accel;
    vector<vec1> m_p;
    std::vector<char*> m_fileNames;
};

#endif // FILELOADER_H
