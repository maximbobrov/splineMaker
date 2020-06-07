#include "globals.h"
#include "tools.h"
#include "fileloader.h"
#include "splineinterpolator.h"

#define NUMCELLS 10

void threadOptimize(int threadIdx, int startIdx, int endIdx);
void display(void);
void init();

FileLoader* fileLoader;

vector<SplineInterpolator*> splineInterpolatorsX;
vector<SplineInterpolator*> splineInterpolatorsY;
vector<SplineInterpolator*> splineInterpolatorsZ;

vector<posVelAccel> currValues;

vector<int> particlesGrid[NUMCELLS][NUMCELLS][NUMCELLS];

void display(void)
{
    double orient_x=0.0;
    double orient_y=0.0;
    double orient_z=5.0;
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glScalef(sc,sc,sc);
    o_x=cos(-ry/20)*cos(-rx/20);
    o_y=cos(-ry/20)*sin(-rx/20);
    o_z=sin(-ry/20);
    orient_x=view_x+o_x;
    orient_y=view_y+o_y;
    orient_z=view_z+o_z;
    gluLookAt(view_x,view_y,view_z,orient_x,orient_y,orient_z,0,0,1);

    size_t numVelVec = 20;
    size_t numSplinePoints = 100;
    size_t drawSplineStep = 100;
    size_t drawPointStep = 1;

    if(!rangeCalculated)
    {
        for (size_t s = 0; s < fileLoader->m_splinesCount; s++)
        {
            if(!splineInterpolatorsX[s]->m_moreMinNumber)
                continue;
            for( size_t i=0; i <= numSplinePoints; i++ )
            {
                double vel = splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints)
                        +    splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints)
                        +    splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints);
                if(minVel > sqrt(vel))
                    minVel = sqrt(vel);
                if(maxVel < sqrt(vel))
                    maxVel = sqrt(vel);

                double accel = splineInterpolatorsX[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsX[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints)
                        +      splineInterpolatorsY[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsY[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints)
                        +      splineInterpolatorsZ[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsZ[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints);
                if(minAccel > accel)
                    minAccel = accel;
                if(maxAccel < accel)
                    maxAccel= accel;
            }
        }
        rangeCalculated = true;
    }

    if(drawPoint)
    {
        /*glPointSize(3);
        glBegin(GL_POINTS);
        for( size_t i = 0; i < currValues.size(); i++ )
        {
            if(currValues[i].isBound)
            {
                double vel = currValues[i].u *  currValues[i].u
                        +    currValues[i].v *  currValues[i].v
                        +    currValues[i].w *  currValues[i].w;
                get_color(sqrt(vel), minVel, maxVel);
                glVertex3f(currValues[i].x, currValues[i].y, currValues[i].z);
            }
        }
        glEnd();

        glPointSize(7);
        glBegin(GL_POINTS);
        glColor3f(1,0,0);
        int num = 0;
        glVertex3f(currValues[num].x, currValues[num].y, currValues[num].z);
        for( size_t i = 0; i <= currValues[num].neighbours.size() - 1; i++ )
        {
            glColor3f(0, 0.2 + i * 0.8 /(currValues[num].neighbours.size()-1),0);
            glVertex3f(currValues[currValues[num].neighbours[i]].x, currValues[currValues[num].neighbours[i]].y, currValues[currValues[num].neighbours[i]].z);
        }
        glEnd();*/
        /*double minP=1e10;
        double maxP=-1e10;
        for (size_t s = 0; s < fileLoader->m_p.size(); s+=100)
        {
            if(!splineInterpolatorsX[s]->m_moreMinNumber)
                continue;
            for( size_t i = 0; i < fileLoader->m_p[s].p.size(); i++ )
            {
                if(currTime == fileLoader->m_p[s].t[i])
                {
                    minP = min(minP, fileLoader->m_p[s].p[i]);
                    maxP = max(maxP, fileLoader->m_p[s].p[i]);
                }
            }
        }*/
        /*glPointSize(3);
        glBegin(GL_POINTS);
        for (size_t s = 0; s < fileLoader->m_p.size(); s+=100)
        {
            if(!splineInterpolatorsX[s]->m_moreMinNumber)
                continue;
            for( size_t i = 0; i < fileLoader->m_p[s].p.size(); i++ )
            {
                //if(currTime == fileLoader->m_p[s].t[i])
                {
                    get_color(fileLoader->m_p[s].p[i], minP, maxP);
                    glVertex3f(fileLoader->m_data[s].x[i], fileLoader->m_data[s].y[i], fileLoader->m_data[s].z[i]);
                }
            }
        }
        glEnd();*/

        glPointSize(5);
        glBegin(GL_POINTS);
        for (size_t s = 0; s < fileLoader->m_vel.size(); s+=100)
        {
            if(!splineInterpolatorsX[s]->m_moreMinNumber)
                continue;
            for( size_t i = 0; i < fileLoader->m_vel[s].x.size(); i++ )
            {
                //if(currTime == fileLoader->m_p[s].t[i])
                {
                    double vel = fileLoader->m_vel[s].x.size() * fileLoader->m_vel[s].x.size()
                            + fileLoader->m_vel[s].y.size() * fileLoader->m_vel[s].y.size()
                            + fileLoader->m_vel[s].z.size() * fileLoader->m_vel[s].z.size();
                    get_color(sqrt(vel), minVel, maxVel);
                    glVertex3f(fileLoader->m_data[s].x[i], fileLoader->m_data[s].y[i], fileLoader->m_data[s].z[i]);
                }
            }
        }
        glEnd();

    }



    for (size_t s = 0; s < fileLoader->m_splinesCount; s+=drawSplineStep)
    {
        if(!splineInterpolatorsX[s]->m_moreMinNumber)
            continue;
        if(drawLineSeg)
        {
            glLineWidth(2);
            glBegin(GL_LINE_STRIP);
            for( size_t i = 0; i < splineInterpolatorsX[s]->m_size; i+=drawPointStep )
            {
                glColor3f(1.0,1.0,1.0);
                glVertex3f(splineInterpolatorsX[s]->m_vec[i], splineInterpolatorsY[s]->m_vec[i], splineInterpolatorsZ[s]->m_vec[i]);
            }
            glEnd();
        }
        if(drawSpline)
        {
            glLineWidth(1);
            glBegin(GL_LINE_STRIP);
            for( size_t i=0; i<=numSplinePoints; i++ )
            {
                double vel = splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints)
                        +    splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints)
                        +    splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints) *  splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints);
                get_color(sqrt(vel), minVel, maxVel);
                glVertex3f(splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numSplinePoints), splineInterpolatorsY[s]->getCoord(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numSplinePoints), splineInterpolatorsZ[s]->getCoord(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numSplinePoints));
            }
            glEnd();
        }
        if(drawVelocity)
        {
            glLineWidth(1);
            glBegin(GL_LINES);
            for( size_t i=0; i<=numVelVec; i++)
            {
                double vel = splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                        +    splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                        +    splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec);
                get_color(sqrt(vel), minVel, maxVel);
                glVertex3f(splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec), splineInterpolatorsY[s]->S_(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numVelVec), splineInterpolatorsZ[s]->S_(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numVelVec));
                glVertex3f(splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)   + scale * 5.0 * splineInterpolatorsX[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                           , splineInterpolatorsY[s]->getCoord(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numVelVec) + scale * 5.0 * splineInterpolatorsY[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                           , splineInterpolatorsZ[s]->getCoord(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numVelVec) + scale * 5.0 * splineInterpolatorsZ[s]->getVel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec));
            }
            glEnd();
        }
        if(drawAcceleration)
        {
            glLineWidth(2);
            glBegin(GL_LINES);
            for( size_t i=0; i<=numVelVec; i++)
            {
                double accel = splineInterpolatorsX[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsX[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                        +      splineInterpolatorsY[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsY[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                        +      splineInterpolatorsZ[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec) *  splineInterpolatorsZ[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec);
                get_color(accel, minAccel, maxAccel);
                glVertex3f( splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec), splineInterpolatorsY[s]->getCoord(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numVelVec), splineInterpolatorsZ[s]->getCoord(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numVelVec));
                glVertex3f( splineInterpolatorsX[s]->getCoord(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)  + scale * splineInterpolatorsX[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                            ,splineInterpolatorsY[s]->getCoord(i*(splineInterpolatorsY[s]->m_size-1)*1.0/numVelVec) + scale * splineInterpolatorsY[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec)
                            ,splineInterpolatorsZ[s]->getCoord(i*(splineInterpolatorsZ[s]->m_size-1)*1.0/numVelVec) + scale * splineInterpolatorsZ[s]->getAccel(i*(splineInterpolatorsX[s]->m_size-1)*1.0/numVelVec));
            }
            glEnd();
        }
    }
    if(drawGridVel)
    {
        for( int k=0; k < min(kNum, NZ); k+=1)
        {
            for( int i=0; i<min(iNum, NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j< min(jNum, NY); j+=1)
                {
                    double vel = fileLoader->m_velField[0][i][j][k] * fileLoader->m_velField[0][i][j][k]
                            +    fileLoader->m_velField[1][i][j][k] * fileLoader->m_velField[1][i][j][k]
                            +    fileLoader->m_velField[2][i][j][k] * fileLoader->m_velField[2][i][j][k];
                    get_color(scale * vel, minVel, maxVel);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * k * 1.0 / NZ);

                    double velip = fileLoader->m_velField[0][i+1][j][k] * fileLoader->m_velField[0][i+1][j][k]
                            +      fileLoader->m_velField[1][i+1][j][k] * fileLoader->m_velField[1][i+1][j][k]
                            +      fileLoader->m_velField[2][i+1][j][k] * fileLoader->m_velField[2][i+1][j][k];
                    get_color(scale * velip, minVel, maxVel);
                    glVertex3f(x0_ + (x1_-x0_) * (i+1) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * k * 1.0 / NZ);
                }
                glEnd();
            }
        }
        for( int i=0; i<min(iNum, NX); i+=1)
        {
            for( int k=0; k < min(kNum, NZ)-1; k+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j<min(jNum, NY); j+=1)
                {
                    double vel = fileLoader->m_velField[0][i][j][k] * fileLoader->m_velField[0][i][j][k]
                            +    fileLoader->m_velField[1][i][j][k] * fileLoader->m_velField[1][i][j][k]
                            +    fileLoader->m_velField[2][i][j][k] * fileLoader->m_velField[2][i][j][k];
                    get_color(scale * vel, minVel, maxVel);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k) * 1.0 / NZ);

                    double velkp = fileLoader->m_velField[0][i][j][k+1] * fileLoader->m_velField[0][i][j][k+1]
                            +      fileLoader->m_velField[1][i][j][k+1] * fileLoader->m_velField[1][i][j][k+1]
                            +      fileLoader->m_velField[2][i][j][k+1] * fileLoader->m_velField[2][i][j][k+1];
                    get_color(scale * velkp, minVel, maxVel);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k+1) * 1.0 / NZ);
                }
                glEnd();
            }
        }
        for( int j=0; j<min(jNum, NY); j+=1)
        {
            for( int i=0; i<min(iNum, NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int k=0; k < min(kNum, NZ); k+=1)
                {
                    double vel = fileLoader->m_velField[0][i][j][k] * fileLoader->m_velField[0][i][j][k]
                            +    fileLoader->m_velField[1][i][j][k] * fileLoader->m_velField[1][i][j][k]
                            +    fileLoader->m_velField[2][i][j][k] * fileLoader->m_velField[2][i][j][k];
                    get_color(scale * vel, minVel, maxVel);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k) * 1.0 / NZ);

                    double velip = fileLoader->m_velField[0][i+1][j][k] * fileLoader->m_velField[0][i+1][j][k]
                            +      fileLoader->m_velField[1][i+1][j][k] * fileLoader->m_velField[1][i+1][j][k]
                            +      fileLoader->m_velField[2][i+1][j][k] * fileLoader->m_velField[2][i+1][j][k];
                    get_color(scale * velip, minVel, maxVel);
                    glVertex3f(x0_ + (x1_-x0_) * (i+1) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k) * 1.0 / NZ);
                }
                glEnd();
            }
        }
    }

    if(drawGridAccel)
    {
        for( int k=0; k < min(kNum, NZ); k+=1)
        {
            for( int i=0; i<min(iNum, NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j< min(jNum, NY); j+=1)
                {
                    double accel = fileLoader->m_accelField[0][i][j][k] * fileLoader->m_accelField[0][i][j][k]
                            +    fileLoader->m_accelField[1][i][j][k] * fileLoader->m_accelField[1][i][j][k]
                            +    fileLoader->m_accelField[2][i][j][k] * fileLoader->m_accelField[2][i][j][k];
                    get_color(scale * accel, minAccel, maxAccel);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * k * 1.0 / NZ);

                    double accelip = fileLoader->m_accelField[0][i+1][j][k] * fileLoader->m_accelField[0][i+1][j][k]
                            +      fileLoader->m_accelField[1][i+1][j][k] * fileLoader->m_accelField[1][i+1][j][k]
                            +      fileLoader->m_accelField[2][i+1][j][k] * fileLoader->m_accelField[2][i+1][j][k];
                    get_color(scale * accelip, minAccel, maxAccel);
                    glVertex3f(x0_ + (x1_-x0_) * (i+1) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * k * 1.0 / NZ);
                }
                glEnd();
            }
        }
        for( int i=0; i<min(iNum, NX); i+=1)
        {
            for( int k=0; k < min(kNum, NZ)-1; k+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j<min(jNum, NY); j+=1)
                {
                    double accel = fileLoader->m_accelField[0][i][j][k] * fileLoader->m_accelField[0][i][j][k]
                            +    fileLoader->m_accelField[1][i][j][k] * fileLoader->m_accelField[1][i][j][k]
                            +    fileLoader->m_accelField[2][i][j][k] * fileLoader->m_accelField[2][i][j][k];
                    get_color(scale * accel, minAccel, maxAccel);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k) * 1.0 / NZ);

                    double accelkp = fileLoader->m_accelField[0][i][j][k+1] * fileLoader->m_accelField[0][i][j][k+1]
                            +      fileLoader->m_accelField[1][i][j][k+1] * fileLoader->m_accelField[1][i][j][k+1]
                            +      fileLoader->m_accelField[2][i][j][k+1] * fileLoader->m_accelField[2][i][j][k+1];
                    get_color(scale * accelkp, minAccel, maxAccel);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k+1) * 1.0 / NZ);
                }
                glEnd();
            }
        }
        for( int j=0; j<min(jNum, NY); j+=1)
        {
            for( int i=0; i<min(iNum, NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int k=0; k < min(kNum, NZ); k+=1)
                {
                    double accel = fileLoader->m_accelField[0][i][j][k] * fileLoader->m_accelField[0][i][j][k]
                            +    fileLoader->m_accelField[1][i][j][k] * fileLoader->m_accelField[1][i][j][k]
                            +    fileLoader->m_accelField[2][i][j][k] * fileLoader->m_accelField[2][i][j][k];
                    get_color(scale * accel, minAccel, maxAccel);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k) * 1.0 / NZ);

                    double accelip = fileLoader->m_accelField[0][i+1][j][k] * fileLoader->m_accelField[0][i+1][j][k]
                            +      fileLoader->m_accelField[1][i+1][j][k] * fileLoader->m_accelField[1][i+1][j][k]
                            +      fileLoader->m_accelField[2][i+1][j][k] * fileLoader->m_accelField[2][i+1][j][k];
                    get_color(scale * accelip, minAccel, maxAccel);
                    glVertex3f(x0_ + (x1_-x0_) * (i+1) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k) * 1.0 / NZ);
                }
                glEnd();
            }
        }
    }

    if(drawGridPressure)
    {
        for( int k=0; k < min(kNum, NZ); k+=1)
        {
            for( int i=0; i<min(iNum, NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j< min(jNum, NY); j+=1)
                {
                    double pressure = fileLoader->m_pressureField[i][j][k];
                    get_color(scale * pressure, minPressure, maxPressure);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * k * 1.0 / NZ);

                    double pressureip = fileLoader->m_pressureField[i+1][j][k];
                    get_color(scale * pressureip, minPressure, maxPressure);
                    glVertex3f(x0_ + (x1_-x0_) * (i+1) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * k * 1.0 / NZ);
                }
                glEnd();
            }
        }
        for( int i=0; i<min(iNum, NX); i+=1)
        {
            for( int k=0; k < min(kNum, NZ)-1; k+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j<min(jNum, NY); j+=1)
                {
                    double pressure =  fileLoader->m_pressureField[i][j][k];
                    get_color(scale * pressure, minPressure, maxPressure);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k) * 1.0 / NZ);

                    double pressurekp =  fileLoader->m_pressureField[i][j][k+1];
                    get_color(scale * pressurekp, minPressure, maxPressure);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k+1) * 1.0 / NZ);
                }
                glEnd();
            }
        }
        for( int j=0; j<min(jNum, NY); j+=1)
        {
            for( int i=0; i<min(iNum, NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int k=0; k < min(kNum, NZ); k+=1)
                {
                    double pressure =  fileLoader->m_pressureField[i][j][k];
                    get_color(scale * pressure, minPressure, maxPressure);
                    glVertex3f(x0_ + (x1_-x0_) * (i) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k) * 1.0 / NZ);

                    double pressureip =  fileLoader->m_pressureField[i+1][j][k];
                    get_color(scale * pressureip, minPressure, maxPressure);
                    glVertex3f(x0_ + (x1_-x0_) * (i+1) * 1.0 / NX,   y0_ + (y1_-y0_) * (j) * 1.0 / NY,   z0_ + (z1_-z0_) * (k) * 1.0 / NZ);
                }
                glEnd();
            }
        }
    }


    glutSwapBuffers();
    if (redr==1) glutPostRedisplay();
}

void m_m(int x,int y)
{
    if (rotate==1)
    {
        rx=rx0+0.1*(x-mx0);
        ry=ry0+0.1*(y-my0);
    }
    glutPostRedisplay();
}

void m_d(int button, int state,int x, int y)
{
    if (state==GLUT_UP)
    {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN)
    {
        rotate=1;
        mx0=x;
        my0=y;
    }
    mouse_x=(1.0*x)/W_WIDTH;
    mouse_y=(W_HEIGHT-(1.0*y))/W_HEIGHT;
    glutPostRedisplay();
}

double loadSettings()
{
    char name[5000];
    char name1[5000];
    FILE *file_data=fopen("settings.txt", "r");
    fscanf(file_data,"%s %s",&name, &name1);
    fclose(file_data);
    return 0;
}

void threadOptimize(int threadIdx, int startIdx, int endIdx)
{
    for (int s = startIdx; s < endIdx; s++)
    {
        if(!splineInterpolatorsX[s]->m_moreMinNumber)
            continue;
        isNoOptimized[threadIdx] += splineInterpolatorsX[s]->optimizeByGrad(1);
        isNoOptimized[threadIdx] += splineInterpolatorsY[s]->optimizeByGrad(1);
        isNoOptimized[threadIdx] += splineInterpolatorsZ[s]->optimizeByGrad(1);
    }
}

void calcEfficient()
{
    for (size_t numToSkip = 0; numToSkip < 4; numToSkip++)
    {
        double averPosError = 0.0;
        double averVelError = 0.0;
        double averAccelError = 0.0;
        double meanVel = 0.0;
        int number = 0;
        for (size_t s = 0; s < fileLoader->m_splinesCount; s++)
        {
            if(!splineInterpolatorsX[s]->m_moreMinNumber)
                continue;
            for( size_t i = numToSkip; i < splineInterpolatorsX[s]->m_size - numToSkip; i++ )
            {
                double posX = splineInterpolatorsX[s]->getCoord(i) - fileLoader->m_data[s].x[i];
                double posY = splineInterpolatorsY[s]->getCoord(i) - fileLoader->m_data[s].y[i];
                double posZ = splineInterpolatorsZ[s]->getCoord(i) - fileLoader->m_data[s].z[i];
                averPosError += posX * posX + posY * posY + posZ * posZ;

                double vecX = splineInterpolatorsX[s]->getVel(i) - fileLoader->m_vel[s].x[i];
                double vecY = splineInterpolatorsY[s]->getVel(i) - fileLoader->m_vel[s].y[i];
                double vecZ = splineInterpolatorsZ[s]->getVel(i) - fileLoader->m_vel[s].z[i];
                averVelError += vecX * vecX + vecY * vecY + vecZ * vecZ;
                meanVel+=sqrt(fileLoader->m_vel[s].x[i]*fileLoader->m_vel[s].x[i] + fileLoader->m_vel[s].y[i]*fileLoader->m_vel[s].y[i] + fileLoader->m_vel[s].z[i]*fileLoader->m_vel[s].z[i]);


                double accelX = splineInterpolatorsX[s]->getAccel(i) - fileLoader->m_accel[s].x[i];
                double accelY = splineInterpolatorsY[s]->getAccel(i) - fileLoader->m_accel[s].y[i];
                double accelZ = splineInterpolatorsZ[s]->getAccel(i) - fileLoader->m_accel[s].z[i];
                averAccelError += accelX * accelX + accelY * accelY + accelZ * accelZ;

                number++;
            }
        }
        printf("skip=%d DiffVel = %f DiffAccel = %f DiffPos = %f \n", numToSkip, sqrt(averVelError/number),sqrt(averAccelError/number),sqrt(averPosError/number));
    }
}

void interpolateInGrid()
{
    printf("Intarpolation started\n");
    for( int xyz=0; xyz<3; xyz+=1)
        for( int i=0; i<NX; i+=1)
            for( int j=0; j<NY; j+=1)
                for( int k=0; k<NZ; k+=1)
                {
                    fileLoader->m_numInCell[i][j][k]=0;
                    fileLoader->m_velField[xyz][i][j][k]=0.0;
                    fileLoader->m_accelField[xyz][i][j][k]=0.0;
                }

    for (size_t i = 0; i < currValues.size(); i++)
    {
        int xIdx = int((currValues[i].x - x0_)/dx);
        int yIdx = int((currValues[i].y - y0_)/dy);
        int zIdx = int((currValues[i].z - z0_)/dz);
        fileLoader->m_velField[0][xIdx][yIdx][zIdx] += currValues[i].u ;
        fileLoader->m_velField[1][xIdx][yIdx][zIdx] += currValues[i].v ;
        fileLoader->m_velField[2][xIdx][yIdx][zIdx] += currValues[i].w ;
        fileLoader->m_accelField[0][xIdx][yIdx][zIdx] += currValues[i].ax ;
        fileLoader->m_accelField[1][xIdx][yIdx][zIdx] += currValues[i].ay ;
        fileLoader->m_accelField[2][xIdx][yIdx][zIdx] += currValues[i].az ;
        fileLoader->m_numInCell[xIdx][yIdx][zIdx]++;
    }
    for( int i=0; i<NX; i+=1)
        for( int j=0; j<NY; j+=1)
            for( int k = 0; k < NZ; k+=1)
                for( int xyz=0; xyz<3; xyz+=1)
                {
                    if(fileLoader->m_numInCell[i][j][k] != 0)
                    {
                        fileLoader->m_velField[xyz][i][j][k] /= fileLoader->m_numInCell[i][j][k];
                        fileLoader->m_accelField[xyz][i][j][k] /= fileLoader->m_numInCell[i][j][k];
                    }
                    fileLoader->m_velFieldCurr[xyz][i][j][k] = fileLoader->m_velField[xyz][i][j][k];
                    fileLoader->m_accelFieldCurr[xyz][i][j][k] = fileLoader->m_accelField[xyz][i][j][k];

                }
    printf("Intarpolation finished\n");
}

void average(int iter)
{
    printf("Average started\n");

    for( size_t tt=0; tt < iter; tt++ )
    {
        for( int xyz=0; xyz<3; xyz+=1)
            for( int i=1; i<NX-1; i+=1)
                for( int j=1; j<NY-1; j+=1)
                    for( int k = 1; k < NZ-1; k+=1)
                    {
                        if(fileLoader->m_velFieldCurr[xyz][i][j][k] > 1e-7)
                            fileLoader->m_velField[xyz][i][j][k] = fileLoader->m_velFieldCurr[xyz][i][j][k];
                        if(fileLoader->m_accelFieldCurr[xyz][i][j][k] > 1e-7)
                            fileLoader->m_accelField[xyz][i][j][k] = fileLoader->m_accelFieldCurr[xyz][i][j][k];
                    }

        for( int xyz=0; xyz<3; xyz+=1)
            for( int i=1; i<NX-1; i+=1)
                for( int j=1; j<NY-1; j+=1)
                    for( int k = 1; k < NZ-1; k+=1)
                    {
                        fileLoader->m_velFieldFiltered[xyz][i][j][k] = (fileLoader->m_velField[xyz][i][j][k]
                                                                        + fileLoader->m_velField[xyz][i+1][j][k]
                                + fileLoader->m_velField[xyz][i-1][j][k]
                                + fileLoader->m_velField[xyz][i][j+1][k]
                                + fileLoader->m_velField[xyz][i][j-1][k]
                                + fileLoader->m_velField[xyz][i][j][k+1]
                                + fileLoader->m_velField[xyz][i][j][k-1])/7.0;
                        fileLoader->m_accelFieldFiltered[xyz][i][j][k] = (fileLoader->m_accelField[xyz][i][j][k]
                                                                          + fileLoader->m_accelField[xyz][i+1][j][k]
                                + fileLoader->m_accelField[xyz][i-1][j][k]
                                + fileLoader->m_accelField[xyz][i][j+1][k]
                                + fileLoader->m_accelField[xyz][i][j-1][k]
                                + fileLoader->m_accelField[xyz][i][j][k+1]
                                + fileLoader->m_accelField[xyz][i][j][k-1])/7.0;
                    }

        for( int xyz=0; xyz<3; xyz+=1)
            for( int i=1; i<NX-1; i+=1)
                for( int j=1; j<NY-1; j+=1)
                    for( int k = 1; k < NZ-1; k+=1)
                    {
                        fileLoader->m_velField[xyz][i][j][k] = (fileLoader->m_velFieldFiltered[xyz][i][j][k]
                                                                + fileLoader->m_velFieldFiltered[xyz][i+1][j][k]
                                + fileLoader->m_velFieldFiltered[xyz][i-1][j][k]
                                + fileLoader->m_velFieldFiltered[xyz][i][j+1][k]
                                + fileLoader->m_velFieldFiltered[xyz][i][j-1][k]
                                + fileLoader->m_velFieldFiltered[xyz][i][j][k+1]
                                + fileLoader->m_velFieldFiltered[xyz][i][j][k-1])/7.0;
                        fileLoader->m_accelField[xyz][i][j][k] = (fileLoader->m_accelFieldFiltered[xyz][i][j][k]
                                                                  + fileLoader->m_accelFieldFiltered[xyz][i+1][j][k]
                                + fileLoader->m_accelFieldFiltered[xyz][i-1][j][k]
                                + fileLoader->m_accelFieldFiltered[xyz][i][j+1][k]
                                + fileLoader->m_accelFieldFiltered[xyz][i][j-1][k]
                                + fileLoader->m_accelFieldFiltered[xyz][i][j][k+1]
                                + fileLoader->m_accelFieldFiltered[xyz][i][j][k-1])/7.0;
                    }
    }
    printf("Average finished\n");
}

void updateCurrValues()
{
    currValues.clear();
    for (size_t s = 0; s < fileLoader->m_splinesCount; s++)
    {
        if(!splineInterpolatorsX[s]->m_moreMinNumber)
            continue;
        for( size_t i = 0; i < splineInterpolatorsX[s]->m_size; i++ )
        {
            if(currTime == splineInterpolatorsX[s]->m_time[i])
                currValues.push_back(posVelAccel(splineInterpolatorsX[s]->getCoord(i), splineInterpolatorsY[s]->getCoord(i), splineInterpolatorsZ[s]->getCoord(i),
                                                 splineInterpolatorsX[s]->getVel(i), splineInterpolatorsY[s]->getVel(i), splineInterpolatorsZ[s]->getVel(i),
                                                 splineInterpolatorsX[s]->getAccel(i), splineInterpolatorsY[s]->getAccel(i), splineInterpolatorsZ[s]->getAccel(i)));
        }
    }
}


void findNeighbors()
{
    double grid_dx = (xmax - xmin) * 1.0 / (NUMCELLS);
    double grid_dy = (ymax - ymin) * 1.0 / (NUMCELLS);
    double grid_dz = (zmax - zmin) * 1.0 / (NUMCELLS);
    double dr = pow(grid_dx*grid_dy*grid_dz,1.0/3)*0.5;
    for (size_t i = 0; i < NUMCELLS; i++)
    {
        for (size_t j = 0; j < NUMCELLS; j++)
        {
            for (size_t k = 0; k < NUMCELLS; k++)
            {
                particlesGrid[i][j][k].clear();
            }
        }
    }
    for (int i = 0; i < currValues.size(); i++)
    {
        int xIdx = int((currValues[i].x - xmin)/grid_dx);
        int yIdx = int((currValues[i].y - ymin)/grid_dy);
        int zIdx = int((currValues[i].z - zmin)/grid_dz);
        particlesGrid[xIdx][yIdx][zIdx].push_back(i);
    }

    //find boundary by x
    for (int i = 0; i < 2 ; i++)
        for (int j = 0; j < NUMCELLS ; j++)
            for (int k = 0; k < NUMCELLS; k++)
            {
                int iIdx = i == 0 ? 0 : NUMCELLS - 1;
                int size = particlesGrid[iIdx][j][k].size();
                if(size == 0)
                    continue;
                double length = pow (((grid_dx * grid_dy * grid_dz) * 1.0 / size), 1.0/ 3.0);
                int targetNumBound = int(size * length / grid_dx);


                for (int ii = 0; ii < size - 1; ii++)
                    for (int jj = 0; jj < size - ii - 1; jj++)
                        if (currValues[particlesGrid[iIdx][j][k].at(jj)].x > currValues[particlesGrid[iIdx][j][k].at(jj+1)].x)
                            swap(currValues[particlesGrid[iIdx][j][k].at(jj)], currValues[particlesGrid[iIdx][j][k].at(jj+1)]);

                if(iIdx == 0)
                    for (int n = 0; n < targetNumBound ; n++)
                        currValues[particlesGrid[iIdx][j][k].at(n)].isBound = true;
                else
                    for (int n = size - 1; n > size - targetNumBound ; n--)
                        currValues[particlesGrid[iIdx][j][k].at(n)].isBound = true;
            }


    /*//find boundary by y
    for (int i = 0; i < NUMCELLS ; i++)
        for (int j = 0; j < 2 ; j++)
            for (int k = 0; k < NUMCELLS; k++)
            {
                int jIdx = j == 0 ? 0 : NUMCELLS - 1;
                int size = particlesGrid[i][jIdx][k].size();
                if(size == 0)
                    continue;
                double length = pow (((grid_dx * grid_dy * grid_dz) * 1.0 / size), 1.0/ 3.0);
                int targetNumBound = int(size * length / grid_dy);


                for (int ii = 0; ii < size - 1; ii++)
                    for (int jj = 0; jj < size - ii - 1; jj++)
                        if (currValues[particlesGrid[i][jIdx][k].at(jj)].y > currValues[particlesGrid[i][jIdx][k].at(jj+1)].y)
                            swap(currValues[particlesGrid[i][jIdx][k].at(jj)], currValues[particlesGrid[i][jIdx][k].at(jj+1)]);

                if(jIdx == 0)
                    for (int n = 0; n < targetNumBound ; n++)
                        currValues[particlesGrid[i][jIdx][k].at(n)].isBound = true;
                else
                    for (int n = size - 1; n > size - targetNumBound ; n--)
                        currValues[particlesGrid[i][jIdx][k].at(n)].isBound = true;
            }*/

    //find boundary by z
    for (int i = 0; i < NUMCELLS ; i++)
        for (int j = 0; j < NUMCELLS ; j++)
            for (int k = 0; k < 2; k++)
            {
                int kIdx = k == 0 ? 0 : NUMCELLS - 1;
                int size = particlesGrid[i][j][kIdx].size();
                if(size == 0)
                    continue;
                double length = pow (((grid_dx * grid_dy * grid_dz) * 1.0 / size), 1.0/ 3.0);
                int targetNumBound = int(size * length / grid_dz);


                for (int ii = 0; ii < size - 1; ii++)
                    for (int jj = 0; jj < size - ii - 1; jj++)
                        if (currValues[particlesGrid[i][j][kIdx].at(jj)].z > currValues[particlesGrid[i][j][kIdx].at(jj+1)].z)
                            swap(currValues[particlesGrid[i][j][kIdx].at(jj)], currValues[particlesGrid[i][j][kIdx].at(jj+1)]);

                if(kIdx == 0)
                    for (int n = 0; n < targetNumBound ; n++)
                        currValues[particlesGrid[i][j][kIdx].at(n)].isBound = true;
                else
                    for (int n = size - 1; n > size - targetNumBound ; n--)
                        currValues[particlesGrid[i][j][kIdx].at(n)].isBound = true;
            }


    for (int s = 0; s < currValues.size(); s++)
    {

        currValues[s].neighbours.clear();
        int xIdx = int((currValues[s].x - xmin)/grid_dx);
        int yIdx = int((currValues[s].y - ymin)/grid_dy);
        int zIdx = int((currValues[s].z - zmin)/grid_dz);
        int im = max(0, xIdx-1);
        int ip = min(NUMCELLS - 1, xIdx+1);
        int jm = max(0, yIdx-1);
        int jp = min(NUMCELLS - 1, yIdx+1);
        int km = max(0, zIdx-1);
        int kp = min(NUMCELLS - 1, zIdx+1);

        for (int i = im; i <= ip ; i++)
            for (int j = jm; j <=jp; j++)
                for (int k = km; k <=kp; k++)
                    for(int n = 0; n < particlesGrid[i][j][k].size(); n++)
                    {
                        if(i == particlesGrid[i][j][k][n])
                            continue;
                        double r2 = (currValues[s].x - currValues[particlesGrid[i][j][k].at(n)].x) * (currValues[s].x - currValues[particlesGrid[i][j][k].at(n)].x)
                                +   (currValues[s].y - currValues[particlesGrid[i][j][k].at(n)].y) * (currValues[s].y - currValues[particlesGrid[i][j][k].at(n)].y)
                                +   (currValues[s].z - currValues[particlesGrid[i][j][k].at(n)].z) * (currValues[s].z - currValues[particlesGrid[i][j][k].at(n)].z);
                        if(r2 < dr * dr)
                        {
                            currValues[s].neighbours.push_back(particlesGrid[i][j][k].at(n));
                        }
                    }

        for (int i = 0; i < currValues[s].neighbours.size() - 1; i++)
            for (int j = 0; j < currValues[s].neighbours.size() - i - 1; j++)
            {
                double r1 = (currValues[s].x - currValues[currValues[s].neighbours.at(j)].x) * (currValues[s].x - currValues[currValues[s].neighbours.at(j)].x)
                        +   (currValues[s].y - currValues[currValues[s].neighbours.at(j)].y) * (currValues[s].y - currValues[currValues[s].neighbours.at(j)].y)
                        +   (currValues[s].z - currValues[currValues[s].neighbours.at(j)].z) * (currValues[s].z - currValues[currValues[s].neighbours.at(j)].z);
                double r2 = (currValues[s].x - currValues[currValues[s].neighbours.at(j+1)].x) * (currValues[s].x - currValues[currValues[s].neighbours.at(j+1)].x)
                        +   (currValues[s].y - currValues[currValues[s].neighbours.at(j+1)].y) * (currValues[s].y - currValues[currValues[s].neighbours.at(j+1)].y)
                        +   (currValues[s].z - currValues[currValues[s].neighbours.at(j+1)].z) * (currValues[s].z - currValues[currValues[s].neighbours.at(j+1)].z);
                if (r1 > r2)
                    swap(currValues[s].neighbours.at(j), currValues[s].neighbours.at(j+1));
            }
    }
}

void saveInFile()
{

    /*for (size_t s = 0; s < fileLoader->m_splinesCount; s++)
    {
        if(!splineInterpolatorsX[s]->m_moreMinNumber)
            continue;
        char filename[64];
        //mkdir("output");
        //sprintf(filename, "output/splines%d_%d.txt", s, splineInterpolatorsX[s]->m_size);

        //mkdir("outputNoisyTracksFreeEnds");
        //mkdir("outputNoisyTracksZeroVelDeriv");
        //mkdir("outputNoisyTracksZeroAccelDeriv");
        mkdir("outputNoisyTracksFreeEnds");
        sprintf(filename, "outputNoisyTracksFreeEnds/splines%d_%d.txt", s, splineInterpolatorsX[s]->m_size);
        FILE *file_data=fopen(filename,"w");
        for( size_t i = 0; i < splineInterpolatorsX[s]->m_size; i++ )
        {
            fprintf(file_data,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", fileLoader->m_data[s].x[i], fileLoader->m_data[s].y[i], fileLoader->m_data[s].z[i]
                    , splineInterpolatorsX[s]->getCoord(i), splineInterpolatorsY[s]->getCoord(i), splineInterpolatorsZ[s]->getCoord(i)
                    , fileLoader->m_vel[s].x[i], fileLoader->m_vel[s].y[i], fileLoader->m_vel[s].z[i]
                    , splineInterpolatorsX[s]->getVel(i), splineInterpolatorsY[s]->getVel(i), splineInterpolatorsZ[s]->getVel(i)
                    , fileLoader->m_accel[s].x[i], fileLoader->m_accel[s].y[i], fileLoader->m_accel[s].z[i]
                    , splineInterpolatorsX[s]->getAccel(i), splineInterpolatorsY[s]->getAccel(i), splineInterpolatorsZ[s]->getAccel(i));
        }
        fclose(file_data);
    }
    printf("Saving done\n");*/

    /*for (size_t s = 0; s < fileLoader->m_splinesCount; s++)
    {
        if(!splineInterpolatorsX[s]->m_moreMinNumber)
            continue;
        char filename[64];
        sprintf(filename, "output/splines%d_%d.txt", s, splineInterpolatorsX[s]->m_size);
        FILE *file_data=fopen(filename,"w");
        for( size_t i = 0; i < splineInterpolatorsX[s]->m_size; i++ )
        {
            fprintf(file_data,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", fileLoader->m_data[s].x[i], fileLoader->m_data[s].y[i], fileLoader->m_data[s].z[i]
                    , splineInterpolatorsX[s]->getCoord(i), splineInterpolatorsY[s]->getCoord(i), splineInterpolatorsZ[s]->getCoord(i)
                    , splineInterpolatorsX[s]->getVel(i), splineInterpolatorsY[s]->getVel(i), splineInterpolatorsZ[s]->getVel(i)
                    , splineInterpolatorsX[s]->getAccel(i), splineInterpolatorsY[s]->getAccel(i), splineInterpolatorsZ[s]->getAccel(i));
        }
        fclose(file_data);
    }
    printf("Saving done\n");*/
    printf("Saving particles in time %d/%d", currTime, fileLoader->m_timeStepCount);
    char filename[64];
    sprintf(filename, "particles%d_%d.txt", currTime, fileLoader->m_timeStepCount);
    FILE *file_data=fopen(filename,"w");
    for (size_t s = 0; s < fileLoader->m_splinesCount; s++)
    {
        for( size_t i = 0; i < splineInterpolatorsX[s]->m_size; i++ )
        {
            if(currTime == splineInterpolatorsX[s]->m_time[i])
            {
                fprintf(file_data,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", fileLoader->m_data[s].x[i], fileLoader->m_data[s].y[i], fileLoader->m_data[s].z[i]
                        , splineInterpolatorsX[s]->getCoord(i), splineInterpolatorsY[s]->getCoord(i), splineInterpolatorsZ[s]->getCoord(i)
                        , splineInterpolatorsX[s]->getVel(i), splineInterpolatorsY[s]->getVel(i), splineInterpolatorsZ[s]->getVel(i)
                        , splineInterpolatorsX[s]->getAccel(i), splineInterpolatorsY[s]->getAccel(i), splineInterpolatorsZ[s]->getAccel(i));
            }
        }
    }
    fclose(file_data);
}



void calcRHS()
{
    int i, j ,k;

    for (k=1;k<NZ-1;k++)
    {
        for (j=1;j<NY-1;j++)
        {
            for (i=1;i<NX-1;i++)
            {
                fileLoader->m_RHS[i][j][k] = -((fileLoader->m_accelField[0][i+1][j][k] - fileLoader->m_accelField[0][i-1][j][k])/(2.0 * dx)
                        + (fileLoader->m_accelField[1][i][j+1][k] - fileLoader->m_accelField[1][i][j-1][k])/(2.0 * dy)
                        + (fileLoader->m_accelField[2][i][j][k+1] - fileLoader->m_accelField[2][i][j][k-1])/(2.0 * dz));
            }
        }
    }
    for (j=1;j<NY-1;j++)
    {
        for (k=1;k<NZ-1;k++)
        {
            fileLoader->m_RHS[0][j][k]=-fileLoader->m_accelField[0][0][j][k];
            fileLoader->m_RHS[NX-1][j][k]=-fileLoader->m_accelField[0][NX-1][j][k];
        }
    }

    for (i=1;i<NX-1;i++)
    {
        for (k=1;k<NZ-1;k++)
        {
            fileLoader->m_RHS[i][0][k]=-fileLoader->m_accelField[1][i][0][k];
            fileLoader->m_RHS[i][NY-1][k]=-fileLoader->m_accelField[1][i][NY-1][k];
        }
    }

    for (i=1;i<NX-1;i++)
    {
        for (j=1;j<NY-1;j++)
        {
            fileLoader->m_RHS[i][j][0]=-fileLoader->m_accelField[2][i][j][0];
            fileLoader->m_RHS[i][j][NZ-1]=-fileLoader->m_accelField[2][i][j][NZ-1];
        }
    }
}

void calcPressure(int in_) //poisson equation at the cell centers
{
    int i,j,k,nn;
    double a,b,c,d,axb_max;
    a=dx*dx*2.0*(1.0/(dx*dx)+ 1.0/(dy*dy)+1.0/(dz*dz));
    b=-1.0*dx*dx/(dx*dx);
    c=-1.0*dx*dx/(dy*dy);
    d=-1.0*dx*dx/(dz*dz);

    axb_max=0.0;
    for (nn=0;nn<in_;nn++)
    {
        calcRHS();
        //neumann bc
        for (j=1;j<NY-1;j++)
        {
            for (k=1;k<NZ-1;k++)
            {
                fileLoader->m_pressureField[0][j][k]=ddx_a(fileLoader->m_pressureField,0,j,k,fileLoader->m_RHS[0][j][k]);//0
                fileLoader->m_pressureField[NX-1][j][k]=ddx_a(fileLoader->m_pressureField,NX-1,j,k,fileLoader->m_RHS[NX-1][j][k]);//0
            }
        }

        for (i=1;i<NX-1;i++)
        {
            for (k=1;k<NZ-1;k++)
            {
                fileLoader->m_pressureField[i][0][k]=ddy_a(fileLoader->m_pressureField,i,0,k,fileLoader->m_RHS[i][0][k]);
                fileLoader->m_pressureField[i][NY-1][k]=ddy_a(fileLoader->m_pressureField,i,NY-1,k,fileLoader->m_RHS[i][NY-1][k]);
            }
        }

        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                fileLoader->m_pressureField[i][j][0]=ddz_a(fileLoader->m_pressureField,i,j,0,fileLoader->m_RHS[i][j][0]);//0
                fileLoader->m_pressureField[i][j][NZ-1]=ddz_a(fileLoader->m_pressureField,i,j,NZ-1,fileLoader->m_RHS[i][j][NZ-1]);//0
            }
        }

        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                for (k=1;k<NZ-1;k++)
                {
                    double d2phi_dx_a,d2phi_dy_a,d2phi_dz_a,d2phi_dx_res,d2phi_dy_res,d2phi_dz_res,phi_a,phi_res;

                    d2phi_dx_a=d2di_a(fileLoader->m_pressureField,i,j,k)/(dx*dx);
                    d2phi_dy_a=d2dj_a(fileLoader->m_pressureField,i,j,k)/(dy*dy);
                    d2phi_dz_a=d2dk_a(fileLoader->m_pressureField,i,j,k)/(dz*dz);

                    d2phi_dx_res=d2di_res(fileLoader->m_pressureField,i,j,k)/(dx*dx);
                    d2phi_dy_res=d2dj_res(fileLoader->m_pressureField,i,j,k)/(dy*dy);
                    d2phi_dz_res=d2dk_res(fileLoader->m_pressureField,i,j,k)/(dz*dz);

                    phi_a=d2phi_dx_a+d2phi_dy_a+d2phi_dz_a;
                    phi_res=d2phi_dx_res+d2phi_dy_res+d2phi_dz_res;
                    fileLoader->m_pressureField[i][j][k]=-(phi_res+fileLoader->m_RHS[i][j][k])/phi_a;
                }
            }
        }

        double pressure_mean=0.0;
        int nnn=0;

        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {

                for (k=1;k<NZ-1;k++)
                {
                    pressure_mean += fileLoader->m_pressureField[i][j][k];
                    nnn++;
                }
            }
        }

        pressure_mean/=nnn;
        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                for (k=1;k<NZ-1;k++)
                {
                    fileLoader->m_pressureField[i][j][k]-=pressure_mean;
                    if(fileLoader->m_pressureField[i][j][k]< minPressure)
                        minPressure = fileLoader->m_pressureField[i][j][k];
                    if(fileLoader->m_pressureField[i][j][k]> maxPressure)
                        maxPressure = fileLoader->m_pressureField[i][j][k];
                }
            }
        }
    }
}

void kb(unsigned char key, int x, int y)
{
    if (key=='1')
    {
        drawPoint = !drawPoint;
    }
    if (key=='2')
    {
        drawLineSeg = !drawLineSeg;
    }
    if (key=='3')
    {
        drawSpline = !drawSpline;
    }
    if (key=='4')
    {
        drawVelocity = !drawVelocity;
    }
    if (key=='5')
    {
        drawAcceleration = !drawAcceleration;
    }
    if (key=='6')
    {
        drawGridVel = !drawGridVel;
        if(drawGridVel)
        {
            drawGridAccel = false;
            drawGridPressure = false;
            interpolateInGrid();
            average(100);
        }
    }
    if (key=='7')
    {
        drawGridAccel = !drawGridAccel;
        if(drawGridAccel)
        {
            drawGridVel = false;
            drawGridPressure = false;
            interpolateInGrid();
            average(100);
        }
    }
    if (key=='8')
    {
        drawGridPressure = !drawGridPressure;
        if(drawGridPressure)
        {
            drawGridVel = false;
            drawGridAccel = false;
            interpolateInGrid();
            average(100);
            calcPressure(10);
        }
    }
    if (key=='0')
    {
        if(currTime < fileLoader->m_timeStepCount)
        {
            currTime+=1;
            updateCurrValues();
            findNeighbors();
        }

        if(drawGridVel || drawGridAccel)
        {
            interpolateInGrid();
            average(100);
        }
        if(drawGridPressure)
        {
            interpolateInGrid();
            average(100);
            calcPressure(10);
        }
        printf("Time %d/%d\n", currTime, fileLoader->m_timeStepCount);
    }
    if (key=='9')
    {
        if(currTime > 1)
        {
            currTime-=1;
            updateCurrValues();
            findNeighbors();
        }
        if(drawGridVel || drawGridAccel)
        {
            interpolateInGrid();
            average(100);
        }
        if(drawGridPressure)
        {
            interpolateInGrid();
            average(100);
            calcPressure(10);
        }
        printf("Time %d/%d\n", currTime, fileLoader->m_timeStepCount);
    }
    if (key=='b')
    {
        average(100);
    }
    if (key=='o')
    {
        iNum++;
        if(iNum>=NX)
            iNum=NX-1;
    }
    if (key=='i')
    {
        iNum--;
        if(iNum<0)
            iNum=0;
    }
    if (key=='l')
    {
        jNum++;
        if(jNum>=NY)
            jNum=NY-1;
    }
    if (key=='k')
    {
        jNum--;
        if(jNum<0)
            jNum=0;
    }

    if (key=='.')
    {
        kNum++;
        if(kNum>=NZ)
            kNum=NZ-1;
    }
    if (key==',')
    {
        kNum--;
        if(kNum<0)
            kNum=0;
    }
    if (key=='[')
    {
        scale/=1.2;
    }
    if (key==']')
    {
        scale*=1.2;
    }
    if (key=='w')
    {
        view_x+=(o_x)*1.5;
        view_y+=(o_y)*1.5;
        view_z+=(o_z)*1.5;
    }
    if (key=='s')
    {
        view_x-=(o_x)*1.5;
        view_y-=(o_y)*1.5;
        view_z-=(o_z)*1.5;
    }
    if (key=='q')
    {
        view_z+=1.1;
    }
    if (key=='e')
    {
        view_z-=1.1;
    }
    if (key=='a')
    {
        double l2=sqrt(o_y*o_y+o_x*o_x);
        view_y+=(o_x)*1.5/l2;
        view_x+=-(o_y)*1.5/l2;
    }
    if (key=='d')
    {
        double    l2=sqrt(o_y*o_y+o_x*o_x);
        view_x+=(o_y)*1.5/l2;
        view_y+=-(o_x)*1.5/l2;
    }

    if(key == 'c')
    {
        for (size_t s = 0; s < fileLoader->m_splinesCount;s++)
        {
            splineInterpolatorsX[s]->optimizeByGrad(1);
            splineInterpolatorsY[s]->optimizeByGrad(1);
            splineInterpolatorsZ[s]->optimizeByGrad(1);
        }
        printf("One minimization step is done\n");
    }
    if(key == 'z')
    {
        calcEfficient();
    }
    if(key == 'x')
    {
        saveInFile();
    }
    if(key == 'v')
    {
        int numNoMin = 1;
        int numToStop = 0;
        double time1=get_time();
        while (numNoMin != 0 && numToStop < 10000)
        {
            numToStop++;
            numNoMin = 0;
            std::vector <std::thread> th_vec;
            th_vec.clear();
            for (int i = 0; i < THREADNUM; ++i)
            {
                int numForOneProc = THREADNUM == 1 ? fileLoader->m_splinesCount : (int)(fileLoader->m_splinesCount / (THREADNUM - 1));
                int startIdx = i * numForOneProc;
                int endIdx =  (i + 1) * numForOneProc;
                if(i == (THREADNUM - 1))
                    endIdx = fileLoader->m_splinesCount;
                isNoOptimized[i] = 0;
                th_vec.push_back(std::thread(threadOptimize,i, startIdx, endIdx));
            }
            for (int i = 0; i < THREADNUM; ++i)
            {
                th_vec.at(i).join();
                numNoMin += isNoOptimized[i];
            }
            printf("Number no minimazed splines = %d\n", numNoMin);
        }
        double time2=get_time();
        printf("Optimization time = %e \n", time2-time1);
    }
    if (key==' ')
    {
        redr=!redr;
    }

    glutPostRedisplay();
}

void init()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective(45.0f, W_WIDTH*1.0/W_HEIGHT, 0.1f, 250.0f);
    glMatrixMode (GL_MODELVIEW);
    loadSettings();
    fileLoader = new FileLoader(PATH);
    fileLoader->findFileInDir();
    fileLoader->loadSplinesFromVlad();//fileLoader->loadSplinesFromDinar();//fileLoader->loadSplinesFromOneFile();//fileLoader->loadSplinesFromNewFormat();//fileLoader->loadSplines();
    printf("%d tracks loaded\n", fileLoader->m_splinesCount);
    splineInterpolatorsX.resize(fileLoader->m_splinesCount);
    splineInterpolatorsY.resize(fileLoader->m_splinesCount);
    splineInterpolatorsZ.resize(fileLoader->m_splinesCount);
    for (size_t i = 0; i < fileLoader->m_splinesCount; i++)
    {
        splineInterpolatorsX[i] = (new SplineInterpolator(fileLoader->m_data[i].x, fileLoader->m_data[i].t));
        splineInterpolatorsY[i] = (new SplineInterpolator(fileLoader->m_data[i].y, fileLoader->m_data[i].t));
        splineInterpolatorsZ[i] = (new SplineInterpolator(fileLoader->m_data[i].z, fileLoader->m_data[i].t));
    }
    printf("Splines calculated\n");

    updateCurrValues();
    //findNeighbors();

    view_x=xmin-0.5*(xmax-xmin);//xmin-0.5*(xmax-xmin);
    view_y=0.5*(ymax-ymin);//ymin-0.5*(ymax-ymin);
    view_z=0.5*(zmax-zmin);//zmin-0.5*(zmax-zmin);
}

int main(int argc, char** argv)
{
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(W_HEIGHT*(1-0)/(1-0),W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();

    delete fileLoader;
    for (size_t i = 0; i < splineInterpolatorsX.size(); i++) {
        delete splineInterpolatorsX[i];
        delete splineInterpolatorsY[i];
        delete splineInterpolatorsZ[i];
    }
}
