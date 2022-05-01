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

void imprimir_vetor(float V[MAX],int n,int p);
// imprime o vetor V de tamanho n; p==0 na posicao horizontal e p==1 na posicao vertical

void imprimir_matriz(float A[MAX][MAX],int l, int c);
// imprime a matriz A de tamanho lxc

void erro_solucao(float A[MAX][MAX],float x[MAX],float d[MAX],int n);
// compara o resultado A.X com D

int main(){
    // INICIALIZAR VARI�VEIS
    bool ciclica;
    int opcao;
    int n;
    float A[MAX][MAX];
    float a[MAX];
    float b[MAX];
    float c[MAX];
    float d[MAX];
    float l[MAX];
    float u[MAX];
    float y[MAX];
    float x[MAX];
    float z[MAX];


    cout<< "-- Bem vindo ao EP1 --"<< endl;
    cout<< "1 - Inserir a matriz tridiagonal" << endl;
    cout<< "2 - Gerar a matriz do enunciado (ciclica)" << endl;
    cout<< "Digite o numero de uma das opcoes acima: ";
    cin>>opcao;

    cout<< "Digite o numero n de linhas e colunas da matriz Anxn: ";
    cin>> n;

    switch (opcao){
        case 1: // Inserir a matriz
                //Entrando com a matriz inteira
            cout<<"Entre com os valores linha a linha: "<< endl;
            for (int i = 1; i <= n; i++){
                cout<<"Entre individualmente com os valores da linha "<<i<<" : "<< endl;
                for (int j = 1; j <= n; j++){
                    cin >> A[i][j];
                }

                //Diagonais a cada linha
                a[i] = A[i][i-1];
                b[i] = A[i][i];
                c[i] = A[i][i+1];
            }

            //Ciclica?
            cout << "A matriz e' ciclica? [S/N] :";
            char cicl;
            cin >> cicl;
            if (cicl == 'S' || cicl == 's'){ciclica = true;}
            else{ciclica = false;}
            
            if (ciclica){
                a[1] = A[1][n];
                c[n] = A[n][1];
            }
            else{
                a[1] = 0;
                c[n] = 0;
                //Para calculo final
                A[1][n] = 0;
                A[n][1] = 0;
            }

            //Entrada de d
            cout << "Digite os valores do termo d: "<<endl;
            cin.clear();
            fflush(stdin);
            for (int i = 1; i <= n; i++){
                cin >> d[i];
            }
            break;

        case 2: // Gerar a matriz
            ciclica = true;
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
            float aux; aux = n;
            a[n] = (2*aux - 1)/(2*aux);
            b[n] = 2;
            c[n] = 1 - a[n];
            d[n] = cos(2*PI);

            // construir a matriz A
            for (int i = 1; i <=n; i++){
                for (int j = 1; j<=n; j++){
                    if(i==j+1){
                        A[i][j]=a[i];
                    } else if(i==j){
                        A[i][j]=b[i];
                    } else if(i+1==j){
                        A[i][j]=c[i];
                    } else if(i==1 && j==n){
                        A[i][j]=a[i];
                    } else if(i==n && j==1){
                        A[i][j]=c[n];
                    }
                }
            }
            break;
        default: cout<< "escolha invalida";
    }

    //Imprimir a matriz

    cout << endl << "A matriz inserida eh a seguinte:" << endl << endl;;
    imprimir_matriz(A,n,n);
    //Mostra das diagonais
    cout << endl << "Os parametros utilizados sao: "<< endl;
    cout << "a[] = " << endl;
    imprimir_vetor(a,n,0);
    cout << "b[] = " << endl;
    imprimir_vetor(b,n,0);
    cout << "c[] = " << endl;
    imprimir_vetor(c,n,0);
    cout << "d[] = " << endl;
    imprimir_vetor(d,n,1);

    // SOLUCAO DO EP

    decomposicaoLU(a,b,c,l,u,n);
    if (ciclica){ // resolver matriz ciclica
        float dn = d[n];
        float aux[MAX];
        float v[MAX];
        for(int i=1;i<=n;i++){
            v[i]=0;
        }
        v[1]=a[1];
        v[n]=c[n-1];

        // Obter y e z solucoes das equacoes Ty=d ; Tz=v
        solucaoLU(y,aux,d,l,u,c,n-1);
        solucaoLU(z,aux,v,l,u,c,n-1);
        // Solucao Xn
        x[n]= (d[n]-c[n]*y[1]-a[n]*y[n-1])/(b[n]-c[n]*z[1]-a[n]*z[n-1]);
        for (int i=n-1;i>0;i--)
            x[i]=y[i]-x[n]*z[i];


    }else // resolver matriz tridiagonal nao ciclica
        solucaoLU(x,y,d,l,u,c,n);

    // Imprimir resposta
    cout << "\nO vetor solucao da matriz A eh:" << endl;
    imprimir_vetor(x,n,1);

    //Calcular resultado e avaliar se bate
    erro_solucao(A,x,d,n);

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

void imprimir_vetor(float V[MAX],int n,int p){
    if (p==0){ // horizontal
        cout << "| ";
        for(int i=1;i<=n;i++){
            printf("%+7.5f ",V[i]);
            if(i!=n)cout << " ";
        };
        cout << " |" << endl;
    } else{// p==1 -> vertical
        for (int i=1;i<=n;i++) {
            cout <<" | ";
            printf("%+7.5f ",V[i]);
            cout <<"| " << endl;
        }
    }
}

void imprimir_matriz(float A[MAX][MAX],int l, int c){
        for (int i=1;i<=l;i++) {
            if(i==l/2)
                 cout << " Matriz A[] = | ";
            else cout << "              | ";
            for (int j=1;j<=c;j++) {
                printf("%+7.5f ",A[i][j]);
        }
        cout << "|" << endl;
    }
}

void erro_solucao(float A[MAX][MAX],float x[MAX],float d[MAX],int n){
    //Calcular resultado e avaliar se bate
    float solucao[MAX];
    float erro[MAX];
    cout<<endl<<endl<<"Comparacao entre solucao encontrada e esperada:"<<endl;
    for (int i = 1; i <= n; i++) {
        solucao[i] = 0;
        for (int j = 1; j <= n; j++) {
            solucao[i] += A[i][j]*x[j];
        }
        if (d[i]>solucao[i]){
            erro[i]=d[i]-solucao[i];
        } else erro[i]=solucao[i]-d[i];
        cout<<"d"<<i<<": encontrada="<<solucao[i]<<"; esperada="<<d[i]<<"; erro="<<erro[i] <<".\n";
    }
}
