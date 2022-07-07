// EP3 MAP3121
// Diogo Vaccaro 8803195
// Otavio Henrique Monteiro 10774159

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>

#define MAX 100
const double PI = 3.141592653589793238463;
char questao = '1'; // Variavel universal para escolha de questao
using namespace std;

double calcula_integral(double a, double b, double T);
double funcao_escolhida(double x, double y);
double solucao_exata(double x, double y);
void solucaoLU(float x[MAX],float y[MAX], float d[MAX],float l[MAX], float u[MAX], float c[MAX], int n);
void decomposicaoLU(float a[MAX],float b[MAX],float c[MAX],float l[MAX],float u[MAX], int n);
void imprimir_vetor(float V[MAX],int n,int p);


int main()
{
    // INICIALIZAR VARIAVEIS
    //  Parametros Gauss (n=2, w=1)
    double T = 1 / sqrt(3); // Abcissa 2 pontos (uma negativa e outra positiva)
    // Matriz Tridiagonal
    int n; // Dimensao da matriz // NOTE: confirmar se este é o mesmo n que é falado no enunciado, parece que sim
    double a[MAX]; // NOTE: confirmar que e' double e nao pode ser float
    double b[MAX];
    double c[MAX];
    double d[MAX];
    double l[MAX];
    double u[MAX];

    cout.precision(17);

    cout << "Bem-Vindo ao EP3 de MAP3121-2022 \n";

    // EXECUCAO
    // Montagem da matriz com produtos internos (integrais)

    

    // Decomposicao LU

    // Solucao da matriz tridiagonal

    // Resultado exato e comparacoes

    // Finalizar
    cout << "\n\nDigite algum caractere para finalizar.\n";
    char end;
    cin >> end;
    return 1;
}

// INTEGRAL SIMPLES (Quadratura Gaussiana de 2 pontos)
double calcula_integral(double a, double b, double T)
{
    double integral = 0;      // valor da integral
    double ba2 = (b - a) / 2; // valor medio do intervalo ab
    double x, y;
    for (int i = 1; i <= 2; i++)
    {
        if (i == 1) // NOTE: tem jeito melhor de colocar o sinal da abcissa mas assim basta
            x = a + ba2 * (-T + 1);
        else
            x = a + ba2 * (T + 1);
        integral += funcao_escolhida(x, 0);
    }
    integral = integral * ba2; // valor final
    return integral;
}

// FUNCAO PARA CADA QUESTAO
double funcao_escolhida(double x, double y)
{
    switch (questao)
    {
    case '1': // f(x) de Validacao 4.2
        return 12 * x * (1 - x) - 2;
    default:
        return 1;
    }
}

double solucao_exata(double x, double y)
{
    switch (questao)
    {
    case '1': // u(x) de Validacao 4.2
        return x * x * (1 - x) * (1 - x);
    default:
        return 1;
    }
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