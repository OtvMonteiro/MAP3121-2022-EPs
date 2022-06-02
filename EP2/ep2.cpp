//EP2 MAP3121
//Diogo Vaccaro 8803195
//Ot�vio Henrique Monteiro 10774159

#include<stdio.h>
#include <iostream>
#include <cmath>
#include<iomanip>

#define MAX 11

using namespace std;

double calcula_integral(int n, double a, double b, double T[MAX],double W[MAX]){
    double integral=0; // valor da integral
    double ba2 = (b-a)/2; // valor medio do intervalo ab
    double x,y;
    for (int i=1;i<=n;i++){ // somatorio peso*funcao(abscissa)
        x = a + ba2*(T[i]+1);
        y = pow(x,3)+1;//funcao_escolhida f(x)=x^3+1
        integral += y*W[i];
    }
    integral = integral * ba2; // valor final
    return integral;
}

int main(){
    // INICIALIZAR VARI�VEIS
    // N� de n�s
    int n;
    int m;
    int i;
    int N;
    double a,b;
    cout<<"Digite o valor de n de pontos (>=1): \n";
    cin>>n;
    cout<<"Digite o limite inferior de integracao a : \n";
    cin>>a;
    cout<<"Digite o limite superior de integracao b : \n";
    cin>>b;

    m = 2*n-1; // maximo grau do polinomio
    double T[MAX]; //vetor de abscissas
    double W[MAX]; //vetor de pesos
    //

    // Preencher T e W com os valores das abscissas e dos pesos para o n determinado (preencher do 1 at� n)
    // para n=2
    if (n==2){
        T[1]=-sqrt(3)/3;
        T[2]=-T[1];
        W[1]=1;
        W[2]=W[1];
    }
    if (n==6){
        T[4]=0.2386191860831969086305017;
        T[5]=0.6612093864662645136613996;
        T[6]=0.9324695142031520278123016;
        T[3]=-T[4];
        T[2]=-T[5];
        T[1]=-T[6];

        W[4]=0.4679139345726910473898703;
        W[5]=0.3607615730481386075698335;
        W[6]=0.1713244923791703450402961;
        W[3]=W[4];
        W[2]=W[5];
        W[1]=W[6];
    }

    // SOLUCAO
    double integral;
    integral = calcula_integral(n,a,b,T,W);



    // Imprimir resposta
    cout << "Para n igual a " << n <<endl;
    cout << "O resultado da Integra��o �:" << integral << endl;

    //Calcular resultado exato


    //Finalizar
    cout << "\n\nDigite algum caractere para finalizar.\n";
    char end;
    cin >> end;
    return 1;
}
