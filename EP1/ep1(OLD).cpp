#include<stdio.h>
#include <iostream>
#include <cmath>

const double PI = 3.141592653589793238463;
#define MAX 100

using namespace std;

// algoritmo para decomposi��o LU de uma matriz tridiagonal A nxn

void solucaoLU(float x[MAX],float y[MAX], float d[MAX],float l[MAX], float u[MAX], float c[MAX], int n);
// Solu��o da equa��o matricial LUx = d

void decomposicaoLU(float a[MAX],float b[MAX],float c[MAX],float l[MAX],float u[MAX], int n);
// decompoe a matriz tridiagonal A nxn em LU

void imprimir_vetor(float V[MAX],int n);
// imprime o vetor V de tamanho n

int main(){
    // INICIALIZAR VARI�VEIS
    bool ciclica = true;
    int respondido = 0; //bool?
    int n;
    //float A[MAX][MAX];
    float a[MAX];
    float b[MAX];
    float c[MAX];
    float d[MAX];
    float l[MAX];
    float u[MAX];
    float y[MAX];
    float x[MAX];
    float z[MAX];

    n=20;
   //Parametros até n-1
    for (int i = 1; i < n; i++)
    {
        float aux = (float) i;
        a[i] = (2*aux - 1)/(4*aux);
        b[i] = 2;
        c[i] = 1 - a[i];
        d[i] = cos((2*PI*i*i)/(n*n));
    }
    //Parametros em i=n
    a[n] = (2*n - 1)/(2*n);
    b[n] = 2;
    c[n] = 1 - a[n];
    d[n] = cos(2*PI);

    

    //Mostra das diagonais
    cout <<endl<< "Os parametros utilizados sao: "<<endl;
    for (int i = 1; i <= n; i++)
    {
        cout<<"a"<<i<<" = "<<a[i]<<endl;
        cout<<"b"<<i<<" = "<<b[i]<<endl;
        cout<<"c"<<i<<" = "<<c[i]<<endl;
        cout<<"d"<<i<<" = "<<d[i]<<endl;

    }
    

    // SOLUCAO DO EP

    decomposicaoLU(a,b,c,l,u,n);
    while (!respondido){
        if (ciclica){ // resolver matriz c�clica
            float dn = d[n];
            float aux[MAX];
            float v[MAX];
            for(int i=1;i<=n;i++){
                v[i]=0;
            }
            v[1]=a[1];
            v[n]=c[n-1];
            respondido = 1;
            
            // Obter y e z solucoes das equacoes Ty=d ; Tz=v
            solucaoLU(y,aux,d,l,u,c,n-1);
            solucaoLU(z,aux,v,l,u,c,n-1);
            // Solucao Xn
            x[n]= (d[n]-c[n]*y[1]-a[n]*y[n-1])/(b[n]-c[n]*z[1]-a[n]*z[n-1]);
            for (int i=n-1;i>0;i--)
                x[i]=y[i]-x[n]*z[i];


        }else{ // resolver matriz normal
            respondido = 1;
            solucaoLU(x,y,d,l,u,c,n);
        }
    }

    // Imprimir resposta
    cout << "O vetor solucao da matriz A eh:\n";
    imprimir_vetor(x,n); //btw o vetor x é vertical, aqui esta' horizontal


    
    //Finalizar
    cout << "\n\nDigite algum caractere para finalizar.\n";
    char end;
    cin >> end;
    return 1;
}

void decomposicaoLU(float a[MAX],float b[MAX],float c[MAX],float l[MAX],float u[MAX],int n){
    u[1]=b[1];
    for(int i=2;i<=n;i++){
        l[i]=a[i]/u[i-1];
        u[i]=b[i]-l[i]*c[i-1];
    }
}

void solucaoLU(float x[MAX],float y[MAX], float d[MAX],float l[MAX], float u[MAX], float c[MAX],int n){
    // Ly = d:
    y[1]=d[1];
    for (int i=2;i<=n;i++)
        y[i]=d[i]-l[i]*y[i-1];
    // Ux = y;
    x[n]=y[n]/u[n];
    for (int i=n-1;i>0;i--)
        x[i]=(y[i]-c[i]*x[i+1])/u[i];
}

void imprimir_vetor(float V[MAX],int n){
    cout << "[ ";
    for(int i=1;i<=n;i++){
        cout<<V[i];
        if(i!=n)cout << " ; ";
    };
    cout << " ]";
}
