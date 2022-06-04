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
   
    //DADOS DE ABCISSAS E PESOS
    double x6[] = {0.0, -0.9324695142031520278123016, -0.6612093864662645136613996, -0.2386191860831969086305017,
                        0.2386191860831969086305017, 0.6612093864662645136613996, 0.9324695142031520278123016}; 
    double w6[] = {0.0, 0.1713244923791703450402961, 0.3607615730481386075698335, 0.4679139345726910473898703,
                        0.4679139345726910473898703, 0.3607615730481386075698335, 0.1713244923791703450402961};
    double x8[] = {0.0, -0.9602898564975362316835609, -0.7966664774136267395915539, -0.5255324099163289858177390, -0.1834346424956498049394761,
                        0.1834346424956498049394761, 0.5255324099163289858177390, 0.7966664774136267395915539, 0.9602898564975362316835609};
    double w8[] = {0.0, 0.1012285362903762591525314, 0.2223810344533744705443560, 0.3137066458778872873379622, 0.3626837833783619829651504,
                        0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314};
    double x10[] = {0.0, -0.9739065285171717200779640, -0.8650633666889845107320967, -0.6794095682990244062343274, -0.4333953941292471907992659, -0.1488743389816312108848260,
                         0.1488743389816312108848260, 0.4333953941292471907992659, 0.6794095682990244062343274, 0.8650633666889845107320967, 0.9739065285171717200779640};
    double w10[] = {0.0, 0.0666713443086881375935688, 0.1494513491505805931457763, 0.2190863625159820439955349, 0.2692667193099963550912269, 0.2955242247147528701738930,
                         0.2955242247147528701738930, 0.2692667193099963550912269, 0.2190863625159820439955349, 0.1494513491505805931457763, 0.0666713443086881375935688};

    // Preencher T e W com os valores das abscissas e dos pesos para o n determinado (preencher do 1 at� n)
    switch (n)
    {
    case 2:
        T[1]=-sqrt(3)/3;
        T[2]=-T[1];
        W[1]=1;
        W[2]=W[1];
        break;
    case 6:
        for (size_t i = 0; i <= n; i++)
        {
           T[i] = x6[i];
           W[i] = w6[i];
        }
        break;
    case 8:
        for (size_t i = 0; i <= n; i++)
        {
           T[i] = x8[i];
           W[i] = w8[i];
        }
        break;
    case 10:
        for (size_t i = 0; i <= n; i++)
        {
           T[i] = x10[i];
           W[i] = w10[i];
        }
        break;
    default:
        cout << "Valor de n inválido";
        char end; cin >> end;
        return 0;
    }

    

    // SOLUCAO
    double integral;
    integral = calcula_integral(n,a,b,T,W);



    // Imprimir resposta
    cout << "Para n igual a " << n <<endl;
    cout << "O resultado da Integracao e':" << integral << endl;

    //Calcular resultado exato


    //Finalizar
    cout << "\n\nDigite algum caractere para finalizar.\n";
    char end;
    cin >> end;
    return 1;
}
