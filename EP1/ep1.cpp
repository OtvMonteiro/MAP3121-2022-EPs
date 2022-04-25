#include<stdio.h>
#include <iostream>
#define MAX 100

using namespace std;

// algoritmo para decomposição LU de uma matriz tridiagonal A nxn

void solucaoLU(float x[MAX],float y[MAX], float d[MAX],float l[MAX], float u[MAX], float c[MAX], int n);
// Solução da equação matricial LUx = d

void decomposicaoLU(float a[MAX],float b[MAX],float c[MAX],float l[MAX],float u[MAX], int n);
// decompoe a matriz tridiagonal A nxn em LU

void imprimir_vetor(float V[MAX],int n);
// imprime o vetor V de tamanho n

int main(){
    // INICIALIZAR VARIÁVEIS
    char ciclica;
    int respondido = 0;
    int n;
    float a[MAX];
    float b[MAX];
    float c[MAX];
    float d[MAX];
    float l[MAX];
    float u[MAX];
    float y[MAX];
    float x[MAX];
    float z[MAX];


    // PREENCHER VETORES a,b,c e d
    cout << "Digite o numero n da matriz tridiagonal A nxn: " ;
    cin >> n;
    a[1]=0; c[n]=0;
    cout << "Em seguida serah pedido para que digite os valores das diagonais da Matriz, termo a termo:\n";
    for(int i=2;i<=n;i++){ // vetor a
        cout << "Digite o valor do termo a" << i << ": ";
        cin >> a[i];
    }
     for(int i=1;i<=n;i++){ // vetor b
        cout << "Digite o valor do termo b" << i << ": ";
        cin >> b[i];
    }
    for(int i=1;i<=n-1;i++){ // vetor c
        cout << "Digite o valor do termo c" << i << ": ";
        cin >> c[i];
    }
     for(int i=1;i<=n;i++){ // vetor d
        cout << "Digite o valor do termo d" << i << ": ";
        cin >> d[i];
    }

    // SOLUCAO DO EP
    decomposicaoLU(a,b,c,l,u,n);
    while (!respondido){
        cout << "A matriz eh ciclica? [S/N] :";
        cin >> ciclica;

        if (ciclica == 'S' || ciclica == 's'){ // resolver matriz cíclica
            float dn = d[n];
            float aux[MAX];
            float v[MAX];
            for(int i=1;i<=n;i++){
                v[i]=0;
            }
            v[1]=a[1];
            v[n]=c[n-1];
            respondido = 1;
            cout << "Digite o valor do termo a1: ";
            cin >> a[1];
            cout << "Digite o valor do termo c" << n << ": ";
            cin >> c[n];
            // Obter y e z solucoes das equacoes Ty=d ; Tz=v
            solucaoLU(y,aux,d,l,u,c,n-1);
            solucaoLU(z,aux,v,l,u,c,n-1);
            // Solucao Xn
            x[n]= (d[n]-c[n]*y[1]-a[n]*y[n-1])/(b[n]-c[n]*z[1]-a[n]*z[n-1]);
            for (int i=n-1;i>0;i--)
                x[i]=y[i]-x[n]*z[i];


        }else if (ciclica == 'N' || ciclica == 'n'){ // resolver matriz normal
            respondido = 1;
            solucaoLU(x,y,d,l,u,c,n);

        }else cout <<"Opcao invalida";
    }

    // Imprimir resposta
    cout << "O vetor solucao da matriz A eh:\n";
    imprimir_vetor(x,n);
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
