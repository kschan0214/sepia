/***********************************************************************
    MGv2.cpp
    written by Dong Zhou
    zhou.dong@gmail.com
    ---------------------
    created:     4.10.13
    last update: 6.15.13
    ------------------------------------------------
    multigrid algorithm for 3D Poisson's equation
        background field removal 
    ------------------------------------------------
    - coarsest grid is has the greatest depth, finest grid is depth 0
    - 
    - 3D image stack is flattened to 1D array
    - matrix size Nx Ny Nz
        in the .bin file, data is organized column-major
    - matrix indexing: x + y*m_size[0] + z*m_size[0]*m_size[1]
    ------------------------------------------------
    - Poisson's equation L u = f 
    - residual equation L e = r, where r = f - L v, e = u - v
    - in the Coarse Grid Correction scheme
        alpha1 and alpha2 are not distinguished, both set to N_GS
        as Gauss-Seidel is used for relaxation
    - 
 
************************************************************************/ 
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include "mex.h"

#define _CORRECTION_ true

#undef _DEBUG_ 

using namespace std;

extern void _main();

class MG {
public:
    MG(void);
    ~MG(void);
    float* fT;          // total magnetic field
    int* mask_;         // mask needs to be fed to the class
    int* interior_;     // internal points, n_peel layers less than the mask

    unsigned int n_vox[10]; // number of voxels 
    float v_size[3];        // voxel size
    int Ncycle;   // number of V cycles at each level, 1 or 2
    float* rho[10];         // source term for the Poisson equation
    float* v[10];           // local field on all grids     
                            //      also error term for residual equation
private:
    int N1, N2;   // number of Gauss-Seidal relaxation iterations
                  // for all depths but the coarsest. pre and post
    int N3;     // extra iteration after all the FMG procedures
    float tol;              // relaxation tolerance
    int depth;     // depth of the V cycle
    unsigned int n_peel;    // throw away the outer n_peel layers

    bool savemode;
    float aux2;         // prefactor 1/(2/dx^2+2/dy^2+2/dz^2)
    float aux3[3];      // 1/dx^2, 1/dy^2, 1/dz^2
    int aux4[10];       // nx * ny 
    int nx[10];         // matrix size for all grids
    int ny[10];
    int nz[10];
    bool* mask;           // boolean mask 
    bool* interior[10];   // interior mask on all grids     
    float* res[10];         // residual on all grids     
                            //      for the V cycle
    unsigned int scale[10];

    void createInterior(bool*,int);
    void initializeRho(void);
    void initializeMem(void);       // v/err, rho and res stack 
    void initializeConstants(void);
    void initializeInterior(void);      // interior mask for coarser grids
    void determineDepth(void);
    void C2F(float* [],int,bool /* = false */);       
    void F2C(float* [],int);
    void updateRes(float*,float*,int);
    void saveData(void);
    float computeDiff(float*, float*, int);
public:
    void setN1(int n){this->N1 = n;};
    void setN2(int n){this->N2 = n;};
    void setN3(int n){this->N3 = n;};
    void setTol(float t){this->tol = t;};
    void setPeel(int n){this->n_peel = n;};
    void setDepth(int d){this->depth = d;};
    void setSavemode(int n){this->savemode = n;};

    void FMG(int);
    void relax(float*,float*, int,int); // Gauss-Seidel
    void V_cycle(float*[],float*, int);
    void loadfT(const char*);
    void loadfL(const char*);
    void loadMask(const char*);
    void readfT(const float*);
    void readMask(const int*);
    void loadParameters(const char*);
    void feedVoxelSize(const float* v){
        for (int i=0;i<3;++i)
            this->v_size[i] = v[i];  };
    void feedMatrixSize(const int*);

    void preProcessing(void);   // prepare the auxiliary variables
    void postProcessing(void);

};
MG::MG(void){
    this->depth = -1;
    this->N1 = 30;        // GS iteration steps
    this->N2 = 100;        
    this->N3 = 100;        
    this->Ncycle = 1;         // V or W cycle
    this->n_peel = 0;
    this->fT = NULL; 
    this->mask_ = NULL; 
    this->mask = NULL; 
    this->interior_ = NULL; 
    this->tol = 0;
    this->savemode = 0;
    this->nx[0]= 0;
    this->ny[0]= 0;
    this->nz[0]= 0;

    int i;
    this->scale[0] = 1;     // scale of finest grid
    for (i=0;i<9;++i)
        this->scale[i+1] = this->scale[i] * 2;
    for (i=0;i<10;++i){
        this->interior[i] = NULL;
        this->rho[i] = NULL;
        this->v[i] = NULL;
        this->res[i] = NULL;
    }
}
MG::~MG(void){
    delete [] this->fT;
    delete [] this->mask_;
    delete [] this->mask;
    delete [] this->interior_;
    for (int i=0;i<=this->depth;++i){
        delete [] this->interior[i];
        delete [] this->rho[i];
        delete [] this->v[i];
        delete [] this->res[i];
    }
}
void MG::loadParameters(const char* filename){

    ifstream fin(filename);

    if (fin){
        cout << "-------------- load parameters ----------------" << endl;
        string tmp;
        while (fin >> tmp){
        //    cout << tmp << endl;
            if (tmp == "matrix_size:" ){
                fin >> this->nx[0];
                fin >> this->ny[0];
                fin >> this->nz[0];
                this->n_vox[0] = this->nx[0]*this->ny[0]*this->nz[0];
            }
            else if (tmp == "voxel_size:"){
                fin >> this->v_size[0];
                fin >> this->v_size[1];
                fin >> this->v_size[2];
            }
            else if (tmp == "fT_name:"){
                fin >> tmp;
                if (tmp.size() > 3)
                    this->loadfT(tmp.c_str());
            }
            else if (tmp == "mask_name:"){
                fin >> tmp;
                if (tmp.size() > 3)
                    this->loadMask(tmp.c_str());
            }
            else if (tmp == "tolerance:")
                fin >> this->tol;
            else if (tmp == "n_peel:")
                fin >> this->n_peel;
            else if (tmp == "depth:")
                fin >> this->depth;
            else if (tmp == "N1:")
                fin >> this->N1;
            else if (tmp == "N2:")
                fin >> this->N2;
            else if (tmp == "N3:")
                fin >> this->N3;
            else if (tmp == "save_mode:")
                fin >> this->savemode;
            else if (tmp == "fL_name:"){
                fin >> tmp;
                if (tmp.size() > 3)
                    this->loadfL(tmp.c_str());
            }
            else
                cout << "unknown parameter" << endl;
        }
        cout << "parameters loaded." << endl;
        fin.close();
    }
    else{
        cout << "error in loading parameters." << endl;
        exit(1);
    }
}
float MG::computeDiff(float* a, float* b, int d){
    float diff = 0;
    for (int i=0;i<this->n_vox[d];++i){
        if (this->interior[d][i] == true){
            diff += abs(a[i]-b[i]);
        }
    }
    return diff;
}
void MG::C2F(float* vv[], int d,bool cor=false){    // d -> d-1
    #ifdef _DEBUG_
        cout << "\tC2F: " << d << ", correction? " << cor << endl;
    #endif
    /* if correction, then add to itself */
    int i,j,k;
    float prefac = 1/ this->scale[d]/this->scale[d];

    for (i=0;i<this->nx[d];++i){    // on coarse grid
    for (j=0;j<this->ny[d];++j){
    for (k=0;k<this->nz[d];++k){
        if (cor){
            vv[d-1][2*i + 2*j*this->nx[d-1] + 2*k*this->aux4[d-1]] 
                    += vv[d][i + j*this->nx[d] + k*this->aux4[d]]*prefac;
        }
        else{
            vv[d-1][2*i + 2*j*this->nx[d-1] + 2*k*this->aux4[d-1]] 
                    =  vv[d][i + j*this->nx[d] + k*this->aux4[d]];
        }
    }
    }
    }
    for (i=0;i<this->nx[d];++i){    // between two coarse grid
    for (j=0;j<this->ny[d];++j){    // points along z 
    for (k=0;k<this->nz[d]-1;++k){
        if (cor){
            vv[d-1][2*i + 2*j*this->nx[d-1] + (2*k+1) *this->aux4[d-1]] 
                    += (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      + vv[d][i + j*this->nx[d] + (k+1)*this->aux4[d]])
                      /2.0 * prefac;
        }
        else{
            vv[d-1][2*i + 2*j*this->nx[d-1] + (2*k+1) *this->aux4[d-1]] 
                    = (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      + vv[d][i + j*this->nx[d] + (k+1)*this->aux4[d]])/2.0;
        }
    }
    }
    }
    for (i=0;i<this->nx[d];++i){      // between two coarse grid
    for (j=0;j<this->ny[d]-1;++j){    // points along y 
    for (k=0;k<this->nz[d];++k){
        if (cor){
            vv[d-1][2*i + (2*j+1)*this->nx[d-1] + 2*k*this->aux4[d-1]] 
                    += (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      + vv[d][i + (j+1)*this->nx[d] + k*this->aux4[d]])
                    / 2.0 * prefac;
        }
        else{
            vv[d-1][2*i + (2*j+1)*this->nx[d-1] + 2*k*this->aux4[d-1]] 
                    = (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      + vv[d][i + (j+1)*this->nx[d] + k*this->aux4[d]])/2.0;
        }
    }
    }
    }
    for (i=0;i<this->nx[d]-1;++i){    // between two coarse grid
    for (j=0;j<this->ny[d];++j){      // points along x
    for (k=0;k<this->nz[d];++k){
        if (cor){
            vv[d-1][(2*i+1) + 2*j*this->nx[d-1] + 2*k*this->aux4[d-1]] 
                    += (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                     + vv[d][(i+1) + j*this->nx[d] + k*this->aux4[d]])
                    / 2.0 * prefac;
        }
        else{
            vv[d-1][(2*i+1) + 2*j*this->nx[d-1] + 2*k*this->aux4[d-1]] 
                    =  (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                     + vv[d][(i+1) + j*this->nx[d] + k*this->aux4[d]])/2.0;
        }
    }
    }
    }
    for (i=0;i<this->nx[d]-1;++i){    // between two coarse grid
    for (j=0;j<this->ny[d]-1;++j){    // points along diagonal 
    for (k=0;k<this->nz[d];++k){      // in x-y plane, assuming uniform
        if (cor){                     // grid size
            vv[d-1][(2*i+1) + (2*j+1)*this->nx[d-1] + 2*k*this->aux4[d-1]] 
                    += (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i+1 + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i + (j+1)*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i+1 + (j+1)*this->nx[d] + k*this->aux4[d]])
                    / 4.0 * prefac;
        }
        else{
            vv[d-1][(2*i+1) + (2*j+1)*this->nx[d-1] + 2*k*this->aux4[d-1]] 
                    =  (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i+1 + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i + (j+1)*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i+1 + (j+1)*this->nx[d] + k*this->aux4[d]])/4.0;
        }
    }
    }
    }
    for (i=0;i<this->nx[d]-1;++i){    // between two coarse grid
    for (j=0;j<this->ny[d];++j){      // points along diagonal 
    for (k=0;k<this->nz[d]-1;++k){    // in x-z plane
        if (cor){
            vv[d-1][2*i+1 + 2*j*this->nx[d-1] + (2*k+1)*this->aux4[d-1]] 
                   += (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i+1 + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][(i+1) + j*this->nx[d] + (k+1)*this->aux4[d]]
                      +vv[d][i + j*this->nx[d] + (k+1)*this->aux4[d]])
                   / 4.0 * prefac;
        }
        else{
            vv[d-1][2*i+1 + 2*j*this->nx[d-1] + (2*k+1)*this->aux4[d-1]] 
                    = (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i+1 + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][(i+1) + j*this->nx[d] + (k+1)*this->aux4[d]]
                      +vv[d][i + j*this->nx[d] + (k+1)*this->aux4[d]])/4.0;
        }
    }
    }
    }
    for (i=0;i<this->nx[d];++i){       // between two coarse grid
    for (j=0;j<this->ny[d]-1;++j){     // points along diagonal 
    for (k=0;k<this->nz[d]-1;++k){     // in y-z plane
        if (cor){
            vv[d-1][2*i + (2*j+1)*this->nx[d-1] + (2*k+1)*this->aux4[d-1]] 
                   += (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i + (j+1)*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i + (j+1)*this->nx[d] + (k+1)*this->aux4[d]]
                      +vv[d][i + j*this->nx[d] + (k+1)*this->aux4[d]])
                   / 4.0 * prefac;
        }
        else{
            vv[d-1][2*i + (2*j+1)*this->nx[d-1] + (2*k+1)*this->aux4[d-1]] 
                   = (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                     +vv[d][i + (j+1)*this->nx[d] + k*this->aux4[d]]
                     +vv[d][i + (j+1)*this->nx[d] + (k+1)*this->aux4[d]]
                     +vv[d][i + j*this->nx[d] + (k+1)*this->aux4[d]])/4.0;
        }
    }
    }
    }
    for (i=0;i<this->nx[d]-1;++i){    // between two coarse grid
    for (j=0;j<this->ny[d]-1;++j){    // points along body diagonal 
    for (k=0;k<this->nz[d]-1;++k){    // 
        if (cor){
            vv[d-1][2*i+1 + (2*j+1)*this->nx[d-1] + (2*k+1)*this->aux4[d-1]] 
                   += (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i+1 + j*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i + (j+1)*this->nx[d] + k*this->aux4[d]]
                      +vv[d][i + j*this->nx[d] + (k+1)*this->aux4[d]]
                      +vv[d][(i+1) + (j+1)*this->nx[d] + k*this->aux4[d]]
                      +vv[d][(i+1) + j*this->nx[d] + (k+1)*this->aux4[d]]
                      +vv[d][(i+1) + (j+1)*this->nx[d] + (k+1)*this->aux4[d]]
                      +vv[d][i + (j+1)*this->nx[d] + (k+1)*this->aux4[d]])
                   / 8.0 * prefac;
        }
        else{
            vv[d-1][2*i+1 + (2*j+1)*this->nx[d-1] + (2*k+1)*this->aux4[d-1]] 
                   = (vv[d][i + j*this->nx[d] + k*this->aux4[d]]
                     +vv[d][i+1 + j*this->nx[d] + k*this->aux4[d]]
                     +vv[d][i + (j+1)*this->nx[d] + k*this->aux4[d]]
                     +vv[d][i + j*this->nx[d] + (k+1)*this->aux4[d]]
                     +vv[d][(i+1) + (j+1)*this->nx[d] + k*this->aux4[d]]
                     +vv[d][(i+1) + j*this->nx[d] + (k+1)*this->aux4[d]]
                     +vv[d][(i+1) + (j+1)*this->nx[d] + (k+1)*this->aux4[d]]
                     +vv[d][i + (j+1)*this->nx[d] + (k+1)*this->aux4[d]])/8.0;
        }
    }
    }
    }
    #ifdef _DEBUG_
    if (!cor){
        FILE *fout;
        fout = fopen("prolong.bin","wb");
        if (fout){
            fwrite(vv[d-1],sizeof(float),this->n_vox[d-1],fout);
            cout << "prolongation saved." << endl;
        }
        fclose(fout);
    }
    #endif
}
void MG::relax(float* vv, float* rr, int d,int N){
    #ifdef _DEBUG_
        cout << "\trelax " << d << ' ' << this->nx[d] << 'x'
                                    << this->ny[d] << 'x'
                                    << this->nz[d] << endl;
    #endif
    int i;
    /* coarsest grid or user choice, iterate until converge */
    if (d == this->depth ){   
        float err = 1;
        float* old = new float [this->n_vox[d]]();
        while (err > this->tol){
            memcpy(old,vv, this->n_vox[d]*sizeof(float));
            for (i=0;i<this->n_vox[d] ;++i)
            if (this->interior[d][i] == true){
                vv[i] = this->aux3[0]* (vv[i+1]+vv[i-1])
                        +this->aux3[1]*(vv[i+this->nx[d]] 
                                +vv[i-this->nx[d]]) 
                        +this->aux3[2]*(vv[i+this->aux4[d]]
                                +vv[i-this->aux4[d]]) 
                        - rr[i];
                vv[i] *= this->aux2;    // prefactor
            }
            err = this->computeDiff(old,this->v[d],d);
//            cout << "update in relax " << err << endl;
        }
        #ifdef _DEBUG_
            cout << "\tconverged." << endl;
        #endif
        delete [] old;
        return;
    }
    /* finer grids, fixed number of iterations */
    for (int k=0;k<N;++k){
        for (i=0;i<this->n_vox[d] ;++i){
        if (this->interior[d][i] == true){
            vv[i] = this->aux3[0]* (vv[i+1]+ vv[i-1])
                    +this->aux3[1]*(vv[i+this->nx[d]] 
                            + vv[i-this->nx[d]]) 
                    +this->aux3[2]*(vv[i+this->aux4[d]]
                            + vv[i-this->aux4[d]]) 
                    - rr[i];
            vv[i] *= this->aux2;
        }
        }
    }
}
void MG::updateRes(float* u,float* f, int d){   // d -> d+1
    #ifdef _DEBUG_
        cout << "\tupdate residual: " << d << endl;
    #endif
    /*   update r = f - L u   */
    for (int i=0;i<this->n_vox[d];++i){
    if (this->interior[d][i] == true){
        this->res[d][i] = f[i] - (
             (u[i+this->aux4[d]]+u[i-this->aux4[d]]-2*u[i])*this->aux3[2]
                +(u[i+this->nx[d]]+u[i-this->nx[d]]-2*u[i])*this->aux3[1]
                +(u[i+1]+ u[i-1]-2*u[i])*this->aux3[0] );
    }
    }
    this->F2C(this->res,d);
}
void MG::FMG(int d){        // d for depth
    #ifdef _DEBUG_
    cout << "FMG, depth: " << d << endl;
    #endif
    if (d == this->depth)   // coarsest grid   
        this->relax(this->v[d],this->rho[d],d,0); //solve fL to convergence
    else{
        this->FMG(d+1);     // go to coarser grid
        #ifdef _DEBUG_
        cout << " still in FMG " << d << endl;
        #endif
        this->C2F(this->v,d+1);   // no correction
        /* note L v = rho is to be solved at depth d */
        for (int i=0;i<this->Ncycle;++i)    // finer grid
            this->V_cycle(this->v,this->rho[d],d);   
    }
}
void MG::V_cycle(float* u[], float* f, int d){
    #ifdef _DEBUG_
        cout << "V cycle " << d << endl;
    #endif
    /*  at level 0, solve L u = f
        at lower levels, solve L e = r   */
    this->relax(u[d],f,d,this->N1);     
    if (d != this->depth){   // not the coarsest grid
        this->updateRes(u[d],f, d);
        fill_n(u[d+1],this->n_vox[d+1],0);    // this initialization
                                                    // might not be 
                                                    // necessary
        this->V_cycle(u,this->res[d+1],d+1);
        this->C2F(u,d+1,_CORRECTION_);   // correct fine grid solution
    }
    this->relax(u[d],f,d,this->N2);
}
void MG::initializeConstants(void){
    int i;
    this->aux2 = 0;     // prefactor 
    for (i=0;i<3;++i){  // used in recursion relation
        this->aux3[i] = 1.0 / this->v_size[i] / this->v_size[i];
        this->aux2 += this->aux3[i];
    }
    this->aux2 = 1/2.0 /this->aux2;
    for (i=0;i<10;++i)  // used in indexing matrix
        this->aux4[i] = this->nx[i] * this->ny[i];
}
void MG::initializeInterior(void){
    #ifdef _DEBUG_
        cout << "init interior" << endl;
    #endif
    int d,i,x,y,z,count;
    bool* tmp = new bool [this->n_vox[1]]();  // tmp mask

    for (d=1;d<=this->depth;++d){    // the finest is already updated
        /* make the coarser mask */
        for (i=0;i<this->n_vox[d];++i){
            x = i % this->nx[d];                // create 3D matrix index 
            y = (i / this->nx[d]) % this->ny[d];        // for the 
            z = i / this->aux4[d];                      // finer grid
            tmp[i] = this->mask[x*this->scale[d]
                                 + y*this->scale[d]*this->nx[0]
                                 + z*this->scale[d]*this->aux4[0]];
        }
        this->createInterior(tmp,d);
    }
    delete [] tmp;
    #ifdef _DEBUG_
        FILE *fout;
        fout = fopen("interior2.bin","wb");
        if (fout){
            fwrite(this->interior[1],sizeof(bool),this->n_vox[1],fout);
            cout << "interior saved." << endl;
        }
        fclose(fout);
    #endif
}
void MG::initializeMem(void){
    #ifdef _DEBUG_
        cout << "init mem" << endl;
    #endif
    for (int d=0;d<=this->depth;++d){
        this->n_vox[d] = this->nx[d]*this->ny[d]*this->nz[d];
        cout << "depth " << d << " voxel number: " 
                        << this->n_vox[d] << endl;
        delete [] this->rho[d];
        this->rho[d] = new float [this->n_vox[d]]();
        delete [] this->v[d];
        this->v[d] = new float [this->n_vox[d]]();
        delete [] this->interior[d];
        this->interior[d] = new bool [this->n_vox[d]]();
        delete [] this->res[d];
        this->res[d] = new float [this->n_vox[d]]();
        // update matrix size for all grids
        this->nx[d+1] = (1+ this->nx[d]) / 2;
        this->ny[d+1] = (1+ this->ny[d]) / 2;
        this->nz[d+1] = (1+ this->nz[d]) / 2;
    }
}
void MG::initializeRho(void){
    #ifdef _DEBUG_
        cout << "init rho" << endl;
    #endif
    /* note scale is set to 1 for all grids 
       rho at each grid is calculated from Laplacian fT
       instead of F2C(rho0)                             */
    int j,x,y,z;
    float* tmpi = new float [this->n_vox[0]]();
    float* tmpj = new float [this->n_vox[0]]();
    float* tmpk = new float [this->n_vox[0]]();
    float tmp;

    for (int d=0;d<=this->depth;++d){
        for (j=0;j<this->n_vox[d];++j){
            if (d == 0)     // for the finest grid
                tmp = this->fT[j];
            else{
                x = j % this->nx[d];                    // create 3D matrix 
                y = (j / this->nx[d]) % this->ny[d];    // index for the 
                z = j / this->aux4[d];                  // finer grid
                tmp = this->fT[x*this->scale[d] 
                               + y*this->scale[d]*this->nx[0]
                               + z*this->scale[d]*this->aux4[0]];
            }
            tmpi[j] = tmp * this->aux3[0];
            tmpj[j] = tmp * this->aux3[1];
            tmpk[j] = tmp * this->aux3[2];
        }

        /* calculate the Laplace */
        for (j=0;j<this->n_vox[d];++j)
        if (this->interior[d][j]){    // Laplace operation
            this->rho[d][j] = tmpk[j+this->aux4[d]]+tmpk[j-this->aux4[d]]
                            +tmpj[j+this->nx[d]]+tmpj[j-this->nx[d]]
                            +tmpi[j+1] + tmpi[j-1]
                            -2*(tmpi[j]+tmpj[j]+tmpk[j]);
        }
    }
    delete [] tmpi;
    delete [] tmpj;
    delete [] tmpk;

    if (this->savemode == 1){
        FILE *fout;
        fout = fopen("rho.bin","wb");
        if (fout){
            fwrite(this->rho[0],sizeof(float),this->n_vox[0],fout);
            cout << "rho saved." << endl;
        }
        fclose(fout);
    }
}
void MG::readfT(const float* pf){
    if (this->n_vox[0] == 0){
        cout << "matrix size undefined." << endl;
        exit(3);
    }
    delete [] this->fT;

    this->fT  = new float [this->n_vox[0]]();
    memcpy(this->fT,pf, this->n_vox[0]*sizeof(float));
}
void MG::readMask(const int* pmask){
    if (this->n_vox[0] == 0){
        cout << "matrix size undefined." << endl;
        exit(3);
    }

    delete [] this->mask_;
    this->mask_ = new int [this->n_vox[0]]();
    memcpy(this->mask_,pmask, this->n_vox[0]*sizeof(int));

    this->mask = new bool [this->n_vox[0]]();
    for (int i=0;i<this->n_vox[0];++i)    // initialize the boolean mask
        this->mask[i] = (this->mask_[i] == 1);
}
void MG::loadMask(const char* f_mask){
    if (this->n_vox[0] == 0){
        cout << "matrix size undefined." << endl;
        exit(3);
    }
    FILE *fin;
    delete [] this->mask_;
    this->mask_ = new int [this->n_vox[0]]();
    fin = fopen(f_mask,"rb");
    if (fin){
        fread(this->mask_,sizeof(int),this->n_vox[0],fin);
        fclose(fin);
        cout << "mask " << f_mask << " loaded." << endl;
        this->mask = new bool [this->n_vox[0]]();
        for (int i=0;i<this->n_vox[0];++i)    // initialize the boolean mask
            this->mask[i] = (this->mask_[i] == 1);
    }
    else{
        cout << "error in read mask." << endl;
        exit(3);
    }
}
void MG::loadfT(const char* f_field){
    if (this->n_vox[0] == 0){
        cout << "matrix size undefined." << endl;
        exit(3);
    }
    FILE *fin;

    delete [] this->fT;
    this->fT  = new float [this->n_vox[0]]();

    fin = fopen(f_field,"rb");
    if (fin){
        fread(this->fT,sizeof(float),this->n_vox[0],fin);
        fclose(fin);
        cout << "total field " << f_field << " loaded." << endl;
    }
    else{
        cout << "error in read fT." << endl;
        exit(3);
    }
}
void MG::loadfL(const char* filename){
    if (this->n_vox[0] == 0){
        cout << "matrix size undefined." << endl;
        exit(3);
    }
    delete [] this->v[0];
    this->v[0] = new float [this->n_vox[0]]();
    FILE *fin;
    fin = fopen(filename,"rb");
    if (fin){
        fread(this->v[0],sizeof(float),this->n_vox[0],fin);
        fclose(fin);
        cout << "local field " << filename << " loaded." << endl;
    }
    else{
        cout << "error in read initial fL." << endl;
        cout << filename << endl;
    }
}
void MG::preProcessing(void){
    if (this->fT == NULL) {
        cout << "field data undefined." << endl;
        exit(3);
    }
    if (this->mask == NULL){
        cout << "mask undefined." << endl;
        exit(3);
    }
    cout << "--------------- pre processing ----------------" << endl;
    if (this->depth < 0)
        this->determineDepth();
    this->initializeMem();
    this->initializeConstants();

    this->createInterior(this->mask,0);
    cout << "peel off " << this->n_peel 
            << " layers from the current mask." << endl;
    for (int i=0;i<this->n_peel;++i){     
        memcpy(this->mask,this->interior[0],sizeof(bool)*this->n_vox[0]);
        this->createInterior(this->mask,0);
    }

    this->initializeInterior();
    this->initializeRho();

    if (this->n_peel != 0){
        stringstream ss;
        ss << this->n_peel;
        string postfix = ss.str() + ".bin";

        FILE *fout;
/*      fout = fopen(("interior_p"+postfix).c_str(),"wb");
        if (fout){
            fwrite(this->interior[0],sizeof(bool),this->n_vox[0],fout);
            cout << "interior mask saved." << endl;
        }
        fclose(fout);
*/
        if (this->mask != NULL){
            int* mask_int = new int [this->n_vox[0]]();
            for (int i=0;i<n_vox[0];++i)   // initialize the integer mask
                mask_int[i] = this->mask[i];
            fout = fopen(("mask_p"+postfix).c_str(),"wb");
            if (fout){
                fwrite(mask_int,sizeof(int),this->n_vox[0], fout);
                cout << "mask saved." << endl;
            }
            fclose(fout);
            delete mask_int;
        }
    }
}
void MG::postProcessing(void){
    // relax on the finest grid
    this->relax(this->v[0],this->rho[0],0,this->N3);
//    this->saveData();
}
void MG::feedMatrixSize(const int m[]){
    this->nx[0] = m[0];
    this->ny[0] = m[1];
    this->nz[0] = m[2];

    int i;
    this->n_vox[0] = 1;
    for (i=0;i<3;++i)
        this->n_vox[0] *= m[i];
}
void MG::determineDepth(void){
    int tmp;
    int scale = 0;
    /* find the smallest matrix dimension */
    if (this->nx[0]<=this->ny[0] && this->nx[0]<=this->nz[0])
        tmp = this->nx[0];
    else if (this->ny[0]<=this->nx[0] && this->ny[0]<=this->nz[0])
        tmp = this->ny[0];
    else
        tmp = this->nz[0];
    if (tmp <= 0 ){
        cout << "Matrix size undefined. " << endl;
        exit(2);
    }
    while (tmp != 1){
        tmp >>= 1;
        ++ scale;
    } 
    cout << scale << endl;
    this->depth = scale>3?scale-3:0;
    cout << "FMG depth = " << this->depth << endl;
}
void MG::F2C(float* vv[], int d){       // d -> d+1
    #ifdef _DEBUG_
        cout << "\tF2C: " << d << endl;
    #endif
    /* take the on-grid points, instead of full weighting */
    int count = 0,x,y,z;

    for (int i=0;i<this->n_vox[d];++i){
        x = i % this->nx[d];                    // create 3D matrix 
        y = (i / this->nx[d]) % this->ny[d];    // index for the 
        z = i / this->aux4[d];                  // finer grid
        if ((x%2==0) && (y%2==0) &&  (z%2==0)){  // if on the coarser grid
            vv[d+1][count] = vv[d][i];
            ++ count ;
        }
    }
    if (count != this->n_vox[d+1])
        cout << "Wrong dimension in F2C: " << count << endl;
}
void MG::createInterior(bool* mask,int d){
    memcpy(this->interior[d],mask,this->n_vox[d]*sizeof(bool));

    int x,y,z;
    for (int i=0;i<this->n_vox[d]; ++i)
    if (mask[i]){
        x = i % this->nx[d];        // create 3D matrix index
        y = (i / this->nx[d]) % this->ny[d];
        z = i / this->aux4[d];
        
        if (x==0||y==0||z==0||x==this->nx[d]-1      // if the point is on 
            ||y==this->ny[d]-1||z==this->nz[d]-1)   // the stack boundary
            this->interior[d][i] = false;
        else if ( (x>0 && mask[i-1]==false) 
                || (x+1<this->nx[d]&& mask[i+1]==false)
                || (y>0 && mask[i-this->nx[d]]==false)
                || (y+1<this->ny[d] 
                        && mask[i+this->nx[d]] ==false)
                || (z>0 && mask[i-this->aux4[d]] ==false)
                || (z+1<this->nz[d]
                        && mask[i+this->aux4[d]] ==false)
                )
            this->interior[d][i] = false;
    }
}
void MG::saveData(void){
    cout << "-------- save results --------" << endl;
    FILE *fout;
    stringstream ss;
    ss << this->n_peel;
    string postfix = ss.str() + ".bin";

    if (this->v[0] != NULL){
        fout = fopen(("fLp"+postfix).c_str(),"wb");
        if (fout){
            fwrite(this->v[0],sizeof(float),this->n_vox[0],fout);
            cout << "fL saved." << endl;
        }
        fclose(fout);
    }
}
    
// ------------------- interface ----------------------------

void mexFunction( int nlhs, mxArray      *plhs[],
		          int nrhs, const mxArray *prhs[]) {
    if (nrhs != 10) {
    mexErrMsgTxt("MEX requires 10 inputs."); }
    
    int N1, N2, N3, peel, depth;
    int *mask, *m_size, n;
    float* fT, *v_size, tol;
    
    tol = (float)(*mxGetPr(prhs[4]));
    depth = *mxGetPr(prhs[5]);
    peel =(int)(*mxGetPr(prhs[6]));
    N1 = (int) (*mxGetPr(prhs[7]));
    N2 = (int) (*mxGetPr(prhs[8]));
    N3 = (int) (*mxGetPr(prhs[9]));

mexPrintf("%g %d %d\n",tol,depth,N1);
mexEvalString("drawnow;");

    int i;
    double *p,*pp;

    m_size = new int [3]();
    v_size = new float [3]();

    p = mxGetPr(prhs[2]);
    pp = mxGetPr(prhs[3]);
    for (i=0;i<3;++i){
        m_size[i] = p[i];
        v_size[i] = pp[i];
    }

mexPrintf("%d %d %d\n",m_size[0],m_size[1],m_size[2]);
mexPrintf("%g %g %g\n",v_size[0],v_size[1],v_size[2]);

    n = mxGetNumberOfElements(prhs[0]);
mexPrintf("%d \n",n);
mexEvalString("drawnow;");

    fT = new float [n]();
    mask = new int [n]();
    p = mxGetPr(prhs[0]);
    pp = mxGetPr(prhs[1]);
    for (i=0;i<n;++i){
        fT[i] = p[i];
        mask[i] = pp[i];
    }

    MG *a = new MG;
    a->feedMatrixSize(m_size);
    a->feedVoxelSize(v_size);
    a->readMask(mask);
    a->readfT(fT);
    a->setDepth(depth);
    a->setN1(N1);
    a->setN2(N2);
    a->setN3(N3);
    a->setTol(tol);
    a->setPeel(peel);


    a->preProcessing();
    a->FMG(0);
    a->postProcessing();

    // associate output
    int dim = a->n_vox[0];

    mxArray *out_m;
    out_m = plhs[0] = mxCreateDoubleMatrix(dim,1,mxREAL);
    double *c;
    c = mxGetPr(out_m);
    for (int i=0;i<dim;++i)
        c[i] = a->v[0][i];

    // free memory
    delete a;
    delete fT;
    delete mask;
    delete m_size;
    delete v_size;

}







