#include "fileloader.h"

FileLoader::vec3::vec3()
{
}

FileLoader::vec3::vec3(double ix, double iy, double iz, int it)
{
    x.push_back(ix);
    y.push_back(iy);
    z.push_back(iz);
    t.push_back(it);
}

FileLoader::vec1::vec1()
{
}

FileLoader::vec1::vec1(double ip,  int it)
{
    p.push_back(ip);
    t.push_back(it);
}

FileLoader::FileLoader(string iPath)
{
    m_path = iPath;
    for (size_t i = 0; i < 10000000; i++)
    {
        m_numbers[i] = -1;
    }
}

int strCompare(const void * a, const void * b)
{
    return strcmp(*(char**)a, *(char**)b);
}

void FileLoader::findFileInDir()
{
    std::string nameFile,nameForOpen;
    std::string nameDir(m_path);
    DIR *mydir = opendir(nameDir.data());
    if(mydir == NULL) {
        perror("opendir");
    }
    struct dirent* entry;
    while ((entry = readdir(mydir)))
    {
        int len = strlen (entry->d_name);
        if (len >= 4) {
            if (strcmp (".dat", entry->d_name + len - 4) == 0)//if (strcmp (".dat", entry->d_name + len - 4) == 0)
            {
                nameFile = entry->d_name;
                nameForOpen = nameDir+ nameFile;
                char *str = new char[nameForOpen.length() + 1];
                strcpy(str, nameForOpen.c_str());
                m_fileNames.push_back(str);
            }
        }
    }
    char* sort_char_array[m_fileNames.size()];
    for (size_t i = 0 ;i < m_fileNames.size();i++)
        sort_char_array[i] = m_fileNames[i];
    qsort (sort_char_array, m_fileNames.size(), sizeof(sizeof(char**)), strCompare);
    for (size_t i =0 ; i < m_fileNames.size();i++)
    {
        m_fileNames[i] = sort_char_array[i];
        printf("FilesName= %s\n",m_fileNames[i]);
    }
    closedir(mydir);
}

void FileLoader::loadSplines()
{
    m_data.clear();
    m_timeStepCount = m_fileNames.size();
    for (size_t i = 0 ; i < m_fileNames.size(); i++)
    {
        FILE *file_data = fopen(m_fileNames[i], "r");
        printf("name=%s\n", m_fileNames[i]);
        char str[128];
        fgets(str, sizeof(str), file_data);
        fgets(str, sizeof(str), file_data);
        fgets(str, sizeof(str), file_data);
        double xx, yy, zz;
        int num;
        while(!feof (file_data))
        {
            fscanf(file_data,"%lf %lf %lf %d", &xx, &yy, &zz, &num);
            if((num >= 0) && (num < 1000000))
            {
                if(m_numbers[num] == -1)
                {
                    m_numbers[num] = m_data.size();
                    m_data.push_back(vec3(xx, yy, zz, i));
                }
                else
                {
                    m_data[m_numbers[num]].x.push_back(xx);
                    m_data[m_numbers[num]].y.push_back(yy);
                    m_data[m_numbers[num]].z.push_back(zz);
                    m_data[m_numbers[num]].t.push_back(i);
                    xmin = xx < xmin ? xx : xmin;
                    xmax = xx > xmax ? xx : xmax;
                    ymin = yy < ymin ? yy : ymin;
                    ymax = yy > ymax ? yy : ymax;
                    zmin = zz < zmin ? zz : zmin;
                    zmax = zz > zmax ? zz : zmax;
                }
            }
        }
        fclose(file_data);
    }
    m_splinesCount = m_data.size();
}

void FileLoader::loadSplinesFromVlad()
{
    m_data.clear();
    m_vel.clear();
    m_accel.clear();
    m_p.clear();
    m_fileNames.resize(10);
    m_timeStepCount = m_fileNames.size();
    for (size_t i = 0 ; i < m_fileNames.size(); i++)
    {
        FILE *file_data = fopen(m_fileNames[i], "r");
        printf("name=%s\n", m_fileNames[i]);
        char str[128];
        //fgets(str, sizeof(str), file_data);
        //fgets(str, sizeof(str), file_data);
        //fgets(str, sizeof(str), file_data);
        double xx, yy, zz;
        double ux, uy, uz;
        double ax, ay, az;
        double pp;
        int num;
        while(!feof (file_data))
        {
            fscanf(file_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &xx, &yy, &zz, &ux, &uy, &uz, &ax, &ay, &az, &pp, &num);
            if((num >= 0) && (num < 1000000) && num % 5 == 0)
            {
                if(m_numbers[num] == -1)
                {
                    m_numbers[num] = m_data.size();
                    m_data.push_back(vec3(xx, yy, zz, i));
                    m_vel.push_back(vec3(ux, uy, uz, i));
                    m_accel.push_back(vec3(ax, ay, az, i));
                    m_p.push_back(vec1(pp,i));
                }
                else
                {
                    m_data[m_numbers[num]].x.push_back(xx);
                    m_data[m_numbers[num]].y.push_back(yy);
                    m_data[m_numbers[num]].z.push_back(zz);
                    m_data[m_numbers[num]].t.push_back(i);

                    m_vel[m_numbers[num]].x.push_back(ux);
                    m_vel[m_numbers[num]].y.push_back(uy);
                    m_vel[m_numbers[num]].z.push_back(uz);
                    m_vel[m_numbers[num]].t.push_back(i);

                    m_accel[m_numbers[num]].x.push_back(ax);
                    m_accel[m_numbers[num]].y.push_back(ay);
                    m_accel[m_numbers[num]].z.push_back(az);
                    m_accel[m_numbers[num]].t.push_back(i);

                    m_p[m_numbers[num]].p.push_back(pp);
                    m_p[m_numbers[num]].t.push_back(i);

                    xmin = xx < xmin ? xx : xmin;
                    xmax = xx > xmax ? xx : xmax;
                    ymin = yy < ymin ? yy : ymin;
                    ymax = yy > ymax ? yy : ymax;
                    zmin = zz < zmin ? zz : zmin;
                    zmax = zz > zmax ? zz : zmax;
                }
            }
        }
        fclose(file_data);
    }
    m_splinesCount = m_data.size();
}

void FileLoader::loadSplinesFromDinar()
{
    m_data.clear();
    m_vel.clear();
    m_accel.clear();
    m_timeStepCount = m_fileNames.size();
    for (size_t i = 0 ; i < m_fileNames.size(); i++)
    {
        FILE *file_data = fopen(m_fileNames[i], "r");
        printf("name=%s\n", m_fileNames[i]);
        char str[128];
        //fgets(str, sizeof(str), file_data);
        //fgets(str, sizeof(str), file_data);
        //fgets(str, sizeof(str), file_data);
        double xx, yy, zz;
        double ux, uy, uz;
        double ax, ay, az;
        int num;
        while(!feof (file_data))
        {
            fscanf(file_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &xx, &yy, &zz, &xx, &yy, &zz, &ux, &uy, &uz, &ax, &ay, &az, &num);
            if((num >= 0) && (num < 1000000))
            {
                if(m_numbers[num] == -1)
                {
                    m_numbers[num] = m_data.size();
                    m_data.push_back(vec3(xx, yy, zz, i));
                    m_vel.push_back(vec3(ux, uy, uz, i));
                    m_accel.push_back(vec3(ax, ay, az, i));
                }
                else
                {
                    m_data[m_numbers[num]].x.push_back(xx);
                    m_data[m_numbers[num]].y.push_back(yy);
                    m_data[m_numbers[num]].z.push_back(zz);
                    m_data[m_numbers[num]].t.push_back(i);

                    m_vel[m_numbers[num]].x.push_back(ux);
                    m_vel[m_numbers[num]].y.push_back(uy);
                    m_vel[m_numbers[num]].z.push_back(uz);
                    m_vel[m_numbers[num]].t.push_back(i);

                    m_accel[m_numbers[num]].x.push_back(ax);
                    m_accel[m_numbers[num]].y.push_back(ay);
                    m_accel[m_numbers[num]].z.push_back(az);
                    m_accel[m_numbers[num]].t.push_back(i);

                    xmin = xx < xmin ? xx : xmin;
                    xmax = xx > xmax ? xx : xmax;
                    ymin = yy < ymin ? yy : ymin;
                    ymax = yy > ymax ? yy : ymax;
                    zmin = zz < zmin ? zz : zmin;
                    zmax = zz > zmax ? zz : zmax;
                }
            }
        }
        fclose(file_data);
    }
    m_splinesCount = m_data.size();
}

void FileLoader::loadSplinesFromNewFormat()
{
    m_data.clear();
    for (size_t i = 0 ; i < m_fileNames.size(); i++)
    {
        m_data.push_back(vec3());
        m_vel.push_back(vec3());
        m_accel.push_back(vec3());
        FILE *file_data = fopen(m_fileNames[i], "r");
        printf("name=%s\n", m_fileNames[i]);
        double xx, yy, zz, xxf, yyf, zzf, vx, vy, vz, ax, ay, az;
        m_numbers[i] = 0;
        while(!feof (file_data))
        {
            int out=fscanf(file_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &xx, &yy, &zz,  &xxf, &yyf, &zzf, &vx, &vy, &vz, &ax, &ay, &az);
            if (out>0)
            {
                m_data[i].x.push_back(xx);
                m_data[i].y.push_back(yy);
                m_data[i].z.push_back(zz);
                m_vel[i].x.push_back(vx);
                m_vel[i].y.push_back(vy);
                m_vel[i].z.push_back(vz);
                m_accel[i].x.push_back(ax);
                m_accel[i].y.push_back(ay);
                m_accel[i].z.push_back(az);
                m_numbers[i]++;
                xmin = xx < xmin ? xx : xmin;
                xmax = xx > xmax ? xx : xmax;
                ymin = yy < ymin ? yy : ymin;
                ymax = yy > ymax ? yy : ymax;
                zmin = zz < zmin ? zz : zmin;
                zmax = zz > zmax ? zz : zmax;
            }
        }
        fclose(file_data);
    }
    m_splinesCount = m_data.size();
}

void FileLoader::loadSplinesFromOneFile()
{
    m_data.clear();
    FILE *file_data = fopen( PATH, "r");
    int lastNum, currNum, readNum;
    currNum = 0;
    lastNum = 0;
    printf("name=%s\n", PATH);

    m_data.push_back(vec3());
    double xx, yy, zz;
    int pointNum;
    int numbers = 0;
    m_timeStepCount = 0;
    while(!feof (file_data))
    {
        lastNum = currNum;
        int out = fscanf(file_data, "%d%*c%d%*c%lf%*c%lf%*c%lf%*c",&currNum, &pointNum, &xx, &yy, &zz);
        if (out>0 && currNum%1 == 0)
        {
            if(pointNum > m_timeStepCount)
                m_timeStepCount = pointNum;
            if(currNum != lastNum)
            {

                m_numbers[numbers] = 0;
                m_data.push_back(vec3(xx, yy, zz, pointNum));
                m_numbers[numbers]++;
                xmin = xx < xmin ? xx : xmin;
                xmax = xx > xmax ? xx : xmax;
                ymin = yy < ymin ? yy : ymin;
                ymax = yy > ymax ? yy : ymax;
                zmin = zz < zmin ? zz : zmin;
                zmax = zz > zmax ? zz : zmax;
                numbers++;
            }
            else
            {
                m_data[numbers].x.push_back(xx);
                m_data[numbers].y.push_back(yy);
                m_data[numbers].z.push_back(zz);
                m_data[numbers].t.push_back(pointNum);
                m_numbers[numbers]++;
                xmin = xx < xmin ? xx : xmin;
                xmax = xx > xmax ? xx : xmax;
                ymin = yy < ymin ? yy : ymin;
                ymax = yy > ymax ? yy : ymax;
                zmin = zz < zmin ? zz : zmin;
                zmax = zz > zmax ? zz : zmax;
            }
        }
    }
    fclose(file_data);
    m_splinesCount = m_data.size();
}
