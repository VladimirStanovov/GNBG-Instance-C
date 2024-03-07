/// Author: Vladimir Stanovov (vladimirstanovov@yandex.ru)
/// Last edited: March 7th, 2024
/// C++ implementation of Generalized Numerical Benchmark Generator (GNBG) instance for GECCO 2024 competition
/// Includes implementation of simple Differential Evolution (DE) with rand/1 strategy and binomial crossover
/// Problem parameters are read from f#.txt files which should be prepared with python script convert.py from f#.mat
/// Competition page: https://competition-hub.github.io/GNBG-Competition/
/// Reference:
/// D. Yazdani, M. N. Omidvar, D. Yazdani, K. Deb, and A. H. Gandomi, "GNBG: A Generalized
///   and Configurable Benchmark Generator for Continuous Numerical Optimization," arXiv prepring	arXiv:2312.07083, 2023.
/// A. H. Gandomi, D. Yazdani, M. N. Omidvar, and K. Deb, "GNBG-Generated Test Suite for Box-Constrained Numerical Global
///   Optimization," arXiv preprint arXiv:2312.07034, 2023.
/// MATLAB and Python versions: https://github.com/Danial-Yazdani/GNBG-Instances
#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <random>

using namespace std;
unsigned globalseed = unsigned(time(NULL));//2024;//
std::mt19937 generator_uni_i(globalseed);
std::mt19937 generator_uni_r(globalseed+100);
std::uniform_int_distribution<int> uni_int(0,32768);
std::uniform_real_distribution<double> uni_real(0.0,1.0);
int IntRandom(int target) {if(target == 0) return 0; return uni_int(generator_uni_i)%target;}
double Random(double minimal, double maximal) {return uni_real(generator_uni_r)*(maximal-minimal)+minimal;}

class GNBG
{
public:
    int FEval;
    int MaxEvals;
    int Dimension;
    int CompNum;
    double MinCoordinate;
    double MaxCoordinate;
    double AcceptanceThreshold;
    double OptimumValue;
    double BestFoundResult;
    double fval;
    double AcceptanceReachPoint;
    double* CompSigma;
    double* Lambda;
    double* OptimumPosition;
    double* FEhistory;
    double* temp;
    double* a;
    double** CompMinPos;
    double** CompH;
    double** Mu;
    double** Omega;
    double*** RotationMatrix;
    GNBG(int func_num);
    ~GNBG();
    double Fitness(double* xvec);
};
GNBG::GNBG(int func_num)
{
    char buffer[15];
    sprintf(buffer,"f%d.txt",func_num);
    ifstream fin(buffer);
    AcceptanceReachPoint = -1;
    FEval = 0;
    fin>>MaxEvals;
    fin>>AcceptanceThreshold;
    fin>>Dimension;
    fin>>CompNum;
    fin>>MinCoordinate;
    fin>>MaxCoordinate;
    FEhistory = new double[MaxEvals];
    a = new double[Dimension];
    temp = new double[Dimension];
    OptimumPosition = new double[Dimension];
    CompSigma = new double[CompNum];
    Lambda = new double[CompNum];
    CompMinPos = new double*[CompNum];
    CompH = new double*[CompNum];
    Mu = new double*[CompNum];
    Omega = new double*[CompNum];
    RotationMatrix = new double**[CompNum];
    for(int i=0;i!=CompNum;i++)
    {
        CompMinPos[i] = new double[Dimension];
        CompH[i] = new double[Dimension];
        Mu[i] = new double[2];
        Omega[i] = new double[4];
        RotationMatrix[i] = new double*[Dimension];
        for(int j=0;j!=Dimension;j++)
            RotationMatrix[i][j] = new double[Dimension];
    }
    for(int i=0;i!=CompNum;i++)
        for(int j=0;j!=Dimension;j++)
            fin>>CompMinPos[i][j];
    for(int i=0;i!=CompNum;i++)
        fin>>CompSigma[i];
    for(int i=0;i!=CompNum;i++)
        for(int j=0;j!=Dimension;j++)
            fin>>CompH[i][j];
    for(int i=0;i!=CompNum;i++)
        for(int j=0;j!=2;j++)
            fin>>Mu[i][j];
    for(int i=0;i!=CompNum;i++)
        for(int j=0;j!=4;j++)
            fin>>Omega[i][j];
    for(int i=0;i!=CompNum;i++)
        fin>>Lambda[i];
    for(int j=0;j!=Dimension;j++)
        for(int k=0;k!=Dimension;k++)
            for(int i=0;i!=CompNum;i++)
                fin>>RotationMatrix[i][j][k];
    fin>>OptimumValue;
    for(int i=0;i!=Dimension;i++)
        fin>>OptimumPosition[i];
}
double GNBG::Fitness(double* xvec)
{
    double res = 0;
    for(int i=0;i!=CompNum;i++)
    {
        for(int j=0;j!=Dimension;j++)
            a[j] = xvec[j] - CompMinPos[i][j];
        for(int j=0;j!=Dimension;j++)
        {
            temp[j] = 0;
            for(int k=0;k!=Dimension;k++)
                temp[j] += RotationMatrix[i][j][k]*a[k]; //matmul rotation matrix and (x - peak position)
        }
        for(int j=0;j!=Dimension;j++)
        {
            if(temp[j] > 0)
                a[j] = exp(log( temp[j])+Mu[i][0]*(sin(Omega[i][0]*log( temp[j]))+sin(Omega[i][1]*log( temp[j]))));
            else if(temp[j] < 0)
                a[j] =-exp(log(-temp[j])+Mu[i][1]*(sin(Omega[i][2]*log(-temp[j]))+sin(Omega[i][3]*log(-temp[j]))));
            else
                a[j] = 0;
        }
        fval = 0;
        for(int j=0;j!=Dimension;j++)
            fval += a[j]*a[j]*CompH[i][j];
        fval = CompSigma[i] + pow(fval,Lambda[i]);
        res = (i == 0)*fval + (i != 0)*min(res,fval);//if first iter then save fval, else take min
    }
    if(FEval > MaxEvals)
        return res;
    FEhistory[FEval] = res;
    BestFoundResult = (FEval == 0)*res + (FEval != 0)*min(res,BestFoundResult);
    if(FEhistory[FEval] - OptimumValue < AcceptanceThreshold && AcceptanceReachPoint == -1)
       AcceptanceReachPoint = FEval;
    FEval++;
    return res;
}
GNBG::~GNBG()
{
    delete a;
    delete temp;
    delete OptimumPosition;
    delete CompSigma;
    delete Lambda;
    for(int i=0;i!=CompNum;i++)
    {
        delete CompMinPos[i];
        delete CompH[i];
        delete Mu[i];
        delete Omega[i];
        for(int j=0;j!=Dimension;j++)
            delete RotationMatrix[i][j];
        delete RotationMatrix[i];
    }
    delete CompMinPos;
    delete CompH;
    delete Mu;
    delete Omega;
    delete RotationMatrix;
    delete FEhistory;
}

class Optimizer
{
public:
    int NInds;
    int NVars;
    int HistorySaveStep;
    int BestHistoryIndex;
    double F;
    double Cr;
    double Left;
    double Right;
    double BestF;
    double tempfit;
    double** Popul;
    double* Fitness;
    double* Trial;
    double* BestHistory;
    Optimizer(int newNInds, GNBG& gnbg);
    ~Optimizer();
    void Run(GNBG& gnbg);
};
Optimizer::Optimizer(int newNInds, GNBG& gnbg)
{
    if(newNInds < 5) return;
    NInds = newNInds;
    NVars = gnbg.Dimension;
    Left = gnbg.MinCoordinate;
    Right = gnbg.MaxCoordinate;
    BestHistoryIndex = 0;
    HistorySaveStep = 1000;
    Popul = new double*[NInds];
    Fitness = new double[NInds];
    Trial = new double[NVars];
    BestHistory = new double[int(double(gnbg.MaxEvals)/double(HistorySaveStep))];
    for(int i=0;i!=NInds;i++)
    {
        Popul[i] = new double[NVars];
        for(int j=0;j!=NVars;j++)
            Popul[i][j] = Random(Left,Right);
    }
}
Optimizer::~Optimizer()
{
    for(int i=0;i!=NInds;i++)
        delete Popul[i];
    delete Popul;
    delete Fitness;
    delete Trial;
    delete BestHistory;
}
void Optimizer::Run(GNBG& gnbg)
{
    for(int i=0;i!=NInds;i++)
    {
        Fitness[i] = gnbg.Fitness(Popul[i]);//if first then save, else take min
        BestF = (i == 0)*Fitness[i] + (i != 0)*min(BestF,Fitness[i]);
        if(gnbg.FEval % HistorySaveStep == 0)
        {
            BestHistory[BestHistoryIndex] = BestF;
            BestHistoryIndex++;
        }
    }
    while(gnbg.FEval < gnbg.MaxEvals)
    {
        for(int i=0;i!=NInds;i++)
        {
            int r1 = IntRandom(NInds-1);
            r1 += (r1 >= i); //so that r1 and i are different
            int r2, r3;
            do
            {
                r2 = IntRandom(NInds);
            } while(r2 == r1 || r2 == i);
            do
            {
                r3 = IntRandom(NInds);
            } while(r3 == r2 || r3 == r1 || r3 == i);
            int jrand = IntRandom(NVars);
            F = 0.5;
            Cr = 0.9;
            for(int j=0;j!=NVars;j++)
            {
                bool cros = (Random(0.0,1.0) < Cr || j == jrand);//if should crossover then perform mutation, else take target vector component
                Trial[j] = cros*(Popul[r1][j] + F*(Popul[r2][j] - Popul[r3][j])) + (!cros)*Popul[i][j];
                Trial[j] = (Trial[j] > Right)*(Popul[i][j]+Right)*0.5 + (Trial[j] <= Right)*Trial[j]; //midpoint-target BCHM
                Trial[j] = (Trial[j] < Left )*(Popul[i][j]+Left )*0.5 + (Trial[j] >= Left )*Trial[j]; //midpoint-target BCHM
            }
            tempfit = gnbg.Fitness(Trial);
            bool repl = (tempfit <= Fitness[i]);//if better then replace, else keep the same
            for(int j=0;j!=NVars;j++)
                Popul[i][j] = repl*Trial[j] + (!repl)*Popul[i][j];
            Fitness[i] = repl*tempfit + (!repl)*Fitness[i];
            BestF = min(BestF,Fitness[i]);
            if(gnbg.FEval % HistorySaveStep == 0)
            {
                BestHistory[BestHistoryIndex] = BestF;
                BestHistoryIndex++;
            }
            if(gnbg.FEval == gnbg.MaxEvals)
                break;
        }
    }
}
int main()
{
    cout.precision(20);
    int NRuns = 31;
    double* Errors = new double[NRuns];
    double* AcceptancePoints = new double[NRuns];
    for(int func_num=1;func_num!=25;func_num++)
    {
        int NNonEmpty = 0;
        double meanAcceptance = 0;
        double meanError = 0;
        double stdAcceptance = 0;
        double stdError = 0;
        for(int run=0;run!=NRuns;run++)
        {
            cout<<"Func: "<<func_num<<"\tRun: "<<run<<endl;
            GNBG gnbg(func_num);
            Optimizer Opt(/*population size*/ 100, gnbg);
            Opt.Run(gnbg);
            if(true) // save for graphs?
            {
                char buffer[100];
                sprintf(buffer,"Res_DE_f%d_r%d.txt",func_num,run);
                ofstream fout_c(buffer);
                fout_c.precision(20);
                for(int i=0;i!=Opt.BestHistoryIndex;i++)
                    fout_c<<Opt.BestHistory[i]-gnbg.OptimumValue<<"\n";
                fout_c.close();
            }
            Errors[run] = gnbg.BestFoundResult - gnbg.OptimumValue;
            AcceptancePoints[run] = gnbg.AcceptanceReachPoint;
            meanError += Errors[run];
            meanAcceptance += AcceptancePoints[run]*(AcceptancePoints[run] != -1);
            NNonEmpty += (AcceptancePoints[run] != -1);
        }
        meanError /= double(NRuns);
        if(NNonEmpty > 0)
            meanAcceptance = meanAcceptance / double(NNonEmpty);
        else
            meanAcceptance = 0;
        for(int run=0;run!=NRuns;run++)
        {
            stdError += (Errors[run]-meanError)*(Errors[run]-meanError);
            stdAcceptance += (AcceptancePoints[run]-meanAcceptance)*(AcceptancePoints[run]-meanAcceptance)*(AcceptancePoints[run] != -1);
        }
        if(NRuns > 1)
            stdError = sqrt(stdError / double(NRuns-1));
        else
            stdError = 0;
        if(NNonEmpty > 1)
            stdAcceptance = sqrt(stdAcceptance / double(NNonEmpty-1));
        else
            stdAcceptance = 0;
        cout<<"Average FE to reach acceptance result: "<<meanAcceptance<<"\t("<<stdAcceptance<<")\n";
        cout<<"Acceptance Ratio: "<<NNonEmpty/NRuns*100<<"\n";
        cout<<"Final result: "<<meanError<<"\t("<<stdError<<")\n";
    }
    delete Errors;
    delete AcceptancePoints;
    return 0;
}
