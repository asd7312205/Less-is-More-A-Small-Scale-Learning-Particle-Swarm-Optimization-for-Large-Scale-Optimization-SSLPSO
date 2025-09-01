#include <iostream>
#include "CEC2013/Header.h"
#include <random>
#include <ctime>
#include <algorithm>
#include <vector>

using namespace std;
Benchmarks* fp=NULL; //选定测试函数
mt19937 rng;//创建随机数
struct Particle{
    double* position;
    double* velocity;
    double fit;
    double distance;
};

bool cmpFit(Particle i,Particle j)
{
    return i.fit<j.fit;
}

bool cmpDistance(Particle i,Particle j)
{
    return i.distance > j.distance;
}

Benchmarks* generateFuncObj(int funcID){
    Benchmarks *fp;
    // run each of specified function in "configure.ini"
    if (funcID==1){
        fp = new F1();
    }else if (funcID==2){
        fp = new F2();
    }else if (funcID==3){
        fp = new F3();
    }else if (funcID==4){
        fp = new F4();
    }else if (funcID==5){
        fp = new F5();
    }else if (funcID==6){
        fp = new F6();
    }else if (funcID==7){
        fp = new F7();
    }else if (funcID==8){
        fp = new F8();
    }else if (funcID==9){
        fp = new F9();
    }else if (funcID==10){
        fp = new F10();
    }else if (funcID==11){
        fp = new F11();
    }else if (funcID==12){
        fp = new F12();
    }else if (funcID==13){
        fp = new F13();
    }else if (funcID==14){
        fp = new F14();
    }else if (funcID==15){
        fp = new F15();
    }else{
        cerr<<"Fail to locate Specified Function Index"<<endl;
        exit(-1);
    }
    return fp;
}



class Population{
public:
    //构造函数
    Population(int particleSize, Benchmarks*fp, int dim, ofstream &write)
    {
        this->fp=fp;
        nowFEs=0;
        allFEs=3000000;
        this->dim =dim;
        this->particleSize=particleSize;
        xMax=fp->getMaxX();//获取测试集函数的约束条件
        xMin=fp->getMinX();
        pe=vector<Particle>(particleSize);
        uniform_real_distribution<double> rd(xMin,xMax);//随机产生位置
        uniform_real_distribution<double> rdv(xMin*0.2,xMax*0.2);//随机产生速度
        for (int i = 0; i < particleSize; ++i) {
            pe[i].position=new double[dim];
            pe[i].velocity=new double[dim];
            for (int j = 0; j < dim; ++j) {
                pe[i].velocity[j]=rdv(rng);
                pe[i].position[j]=rd(rng);
            }
            pe[i].fit=fp->compute(pe[i].position);
        }
        nowFEs+=particleSize;
        sort(pe.begin(),pe.end(), cmpFit);
        write<<nowFEs<<" "<<pe[0].fit<<endl;
        pe[0].distance=0;
        closeIndex=1;
        for (int i = 1; i < particleSize; ++i) {
            pe[i].distance=0;
            for (int j = 0; j < dim; ++j) {
                pe[i].distance+=pow((pe[i].position[j]-pe[0].position[j]),2);
            }
            pe[i].distance= sqrt(pe[i].distance);
            if(pe[i].distance < pe[closeIndex].distance)
            {
                closeIndex=i;
            }
        }

        flags=new bool[6];
        for (int i = 0; i < 6; ++i) {
            flags[i]= true;
        }
    }
    //析构函数
    ~Population()
    {
        for (int i = 0; i < particleSize; ++i) {
            delete [] pe[i].position;
            delete [] pe[i].velocity;
        }
        delete [] flags;

    }

    double finalResult()
    {
        return pe[0].fit;
    }

    void recordKeyPoint(ofstream &write)
    {
        for (int i = 0; i < 6; ++i) {
            if(nowFEs > (i+1)*(allFEs/6) )
            {
                if(flags[i])
                {
                    write<<nowFEs<<" "<<pe[0].fit<<endl;
                    flags[i]= false;
                    break;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                break;
            }
        }
    }

    void SSLPSO(ofstream &write)
    {
        uniform_real_distribution<double> rdr(0,1);
        int e1,e2,sw,worstIndex;
        int xp1,xp2;
        bool diversityFlag;
        double fai;

        while(nowFEs < allFEs)
        {
            diversityFlag= true;
            fai=0.4-0.2* ((double)nowFEs/(double)allFEs);
            particleSize= ceil(1000-800* sqrt((double)nowFEs/(double)allFEs));
            worstIndex=particleSize-1;
            xp1=0.1*particleSize;
            xp2=0.6*particleSize;

            for (int j = 0; j < dim; ++j) {
                e1 = rand()%xp1;

                do {
                    e2 = rand()%xp2;
                } while (e2 == e1);
                if (e1 > e2) {
                    sw = e1;
                    e1 = e2;
                    e2 = sw;
                }
                pe[worstIndex].velocity[j] = rdr(rng) * pe[worstIndex].velocity[j] +
                                             rdr(rng) * (pe[e1].position[j] - pe[worstIndex].position[j]) +
                                             fai * rdr(rng) * (pe[e2].position[j] - pe[worstIndex].position[j]);
                pe[worstIndex].position[j] = pe[worstIndex].position[j] + pe[worstIndex].velocity[j];
                if (pe[worstIndex].position[j] > xMax) {
                    pe[worstIndex].position[j] = xMax;
                } else if (pe[worstIndex].position[j] < xMin) {
                    pe[worstIndex].position[j] = xMin;
                }
            }

            pe[worstIndex].fit=fp->compute(pe[worstIndex].position);
            nowFEs++;

            if(pe[worstIndex].fit < pe[0].fit)
            {
                for (int i = 0; i < particleSize; ++i) {
                    pe[i].distance=0;
                    for (int j = 0; j < dim; ++j) {
                        pe[i].distance+=pow((pe[i].position[j]-pe[worstIndex].position[j]),2);
                    }
                    pe[i].distance= sqrt(pe[i].distance);
                }

            }
            else {
                pe[worstIndex].distance=0;
                for (int j = 0; j < dim; ++j) {
                    pe[worstIndex].distance += pow((pe[worstIndex].position[j] - pe[0].position[j]), 2);
                }
                pe[worstIndex].distance = sqrt(pe[worstIndex].distance);
            }

            sort(pe.begin(),pe.begin()+particleSize, cmpDistance);
            closeIndex=particleSize-2;
            for (int i = particleSize-2; i >particleSize/2 ; i--) {
                if(pe[i].fit > pe[particleSize/2].fit)
                {
                    diversityFlag= false;
                    break;
                }
            }

            if(diversityFlag)
            {
                for (int j = 0; j < dim; ++j) {
                    e1 = rand()%xp1;
                    do {
                        e2 = rand()%xp2;
                    } while ( e2 == e1 );
                    if (e1 > e2) {
                        sw = e1;
                        e1 = e2;
                        e2 = sw;
                    }

                    pe[closeIndex].velocity[j] = rdr(rng) * pe[closeIndex].velocity[j] +
                                                 rdr(rng) * (pe[e1].position[j] - pe[closeIndex].position[j]) +
                                                 fai * rdr(rng) * (pe[e2].position[j] - pe[closeIndex].position[j]);
                    pe[closeIndex].position[j] = pe[closeIndex].position[j] + pe[closeIndex].velocity[j];
                    if (pe[closeIndex].position[j] > xMax) {
                        pe[closeIndex].position[j] = xMax;
                    } else if (pe[closeIndex].position[j] < xMin) {
                        pe[closeIndex].position[j] = xMin;
                    }


                }
                pe[closeIndex].fit=fp->compute(pe[closeIndex].position);
                nowFEs++;

                if(pe[closeIndex].fit < pe[particleSize-1].fit)
                {
                    for (int i = 0; i < particleSize; ++i) {
                        pe[i].distance=0;
                        for (int j = 0; j < dim; ++j) {
                            pe[i].distance+=pow((pe[i].position[j]-pe[closeIndex].position[j]),2);
                        }
                        pe[i].distance= sqrt(pe[i].distance);

                    }
                }
                else {
                    pe[closeIndex].distance=0;
                    for (int j = 0; j < dim; ++j) {
                        pe[closeIndex].distance += pow((pe[closeIndex].position[j] - pe[particleSize-1].position[j]), 2);
                    }
                    pe[closeIndex].distance = sqrt(pe[closeIndex].distance);
                }

            }

            sort(pe.begin(),pe.begin()+particleSize, cmpFit);
            recordKeyPoint(write);
        }

    }
private:
    vector<Particle> pe;
    int dim;
    double xMax;
    double xMin;
    int bestIndex;
    int particleSize;
    int nowFEs;
    int allFEs;
    bool* flags;
    Benchmarks*fp;
    int closeIndex;

};



int main() {
    srand((unsigned )time(NULL));
    int dim = 1000;//维度
    unsigned funToRun[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    unsigned funNum =15;
    ofstream write;
    int particleSize=1000;

    for (int cnum = 0;cnum < funNum; ++cnum) {
        if (12 == cnum || 13 == cnum) {
            dim = 905;
        } else {
            dim = 1000;
        }

        string name1 = "function";
        string name2;
        name2= to_string(cnum+1);
        name1+=name2;
        name2=".txt";
        name1+=name2;
        write.open(name1);

        fp = generateFuncObj(funToRun[cnum]);
        fp->setDimension(dim);
        rng.seed(random_device()());//创建随机数
        Population pn(particleSize, fp, dim, write);
        pn.SSLPSO(write);
        write << endl;
        write << "function" << cnum + 1 << ": " << pn.finalResult() << endl;
        write.close();
    }
    delete fp;
    return 0;
}
