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
double funcao_escolhida(double x);
double solucao_exata(double x, double y);
void solucaoLU(double x[MAX], double y[MAX], double d[MAX], double l[MAX], double u[MAX], double c[MAX], int n);
void decomposicaoLU(double a[MAX], double b[MAX], double c[MAX], double l[MAX], double u[MAX], int n);
void imprimir_vetor(double V[MAX], int n, int p);

int main()
{
    // INICIALIZAR VARIAVEIS
    //  Parametros Gauss (n=2, w=1)
    double T = 1 / sqrt(3); // Abcissa 2 pontos (uma negativa e outra positiva)
    // Matriz Tridiagonal
    int n;         // Dimensao da matriz // NOTE: confirmar se este é o mesmo n que é falado no enunciado, parece que sim
    double a[MAX]; // NOTE: confirmar que e' double e nao pode ser float
    double b[MAX];
    double c[MAX];
    double d[MAX];
    double l[MAX];
    double u[MAX];
    double x[MAX];
    double y[MAX];

    cout.precision(17);

    cout << "Bem-Vindo ao EP3 de MAP3121-2022 \n";
    // NOTE: Acredito que querem que façamos um loop com os diferentes valores, inicialmente assim ta bom
    cout << "Digite o numero n de linhas e colunas da matriz Anxn: ";
    cin >> n;

    // EXECUCAO
    // Montagem da matriz com produtos internos (integrais)
    double h = 1 / ((double)n + 1);
    for (int i = 1; i <= n; i++) // montado a cada linha (inclui poucos pontos desnecessarios, com lixo)
    {
        double aux = (double)i;
        double xi = aux * h;
        double xi_ant = (aux - 1) * h;
        double xi_prox = (aux + 1) * h;
        b[i] = 2 / h;
        c[i] = -1 / h;
        a[i + 1] = c[i];
        d[i] = calcula_integral(xi_ant, xi_prox, T); // TODO: corrigir erro, provavelmente aqui
    }
    imprimir_vetor(a, n + 1, 1);
    cout << "-------" << endl;
    imprimir_vetor(b, n, 1);
    cout << "-------" << endl;
    imprimir_vetor(c, n, 1);
    cout << "-------" << endl;
    imprimir_vetor(d, n, 1);
    cout << "-------" << endl;

    // Decomposicao LU
    decomposicaoLU(a, b, c, l, u, n);

    // Solucao da matriz tridiagonal
    solucaoLU(x, y, d, l, u, c, n);
    // vetor x é o vetor de alfas
    cout << "-------" << endl;
    imprimir_vetor(x, n, 1);
    cout << "-------" << endl;

    // Resultado exato e comparacoes
    for (int i = 1; i <= n; i++) // montado a cada linha (inclui poucos pontos desnecessarios, com lixo)
    {
        double u_barra = 0.0;
        double aux_i = (double)i;
        double xi = aux_i * h;
        for (int j = 1; j <= n; j++)
        {
            // TODO: refatorar, esta no momento so' pegando
            double aux_j = (double)j;
            double xj_ant = (aux_j - 1) * h;
            double phi_j = (xi - xj_ant) / h;
            u_barra += phi_j * x[j];
        }

        double u_exato = solucao_exata(xi, 0);

        cout << "u" << i << ": encontrado=" << u_barra << "; exato=" << u_exato << "\n";
    }

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
    for (int i = 1; i <= 2; i++) // 2 Pontos
    {
        if (i == 1)
        { // NOTE: tem jeito melhor de colocar o sinal da abcissa mas assim basta
            x = a + ba2 * (-T + 1);
        }
        else
        {
            x = a + ba2 * (T + 1);
        }

        double phi = (x - a) / ba2; //(x-x{i-1} /h) em [x_{i-1}, x_{i}] // NOTE: fala pra integrar em cada subintervalo de nos consecutivos, eh so isso? ou o outro intervalo tamb precisa ser considerado?
        integral += phi * funcao_escolhida(x);
        // TODO: considerar casos com k!=1 (talvez so' multiplicar aqui baste)
        // TODO: avaliar necessidade e como implementar q(x)
    }
    integral = integral * ba2; // valor final
    return integral;
}

// FUNCAO PARA CADA QUESTAO
double funcao_escolhida(double x)
{
    switch (questao)
    {
    case '1': // f(x) de Validacao 4.2
        return (12 * x * (1 - x)) - 2;
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

// double phi_i(double x, double x_ant, double h)
// {
//     return (x - x_ant) / h; // em [x_{i-1}, x_{i}] // NOTE: fala pra integrar em cada subintervalo de nos consecutivos, eh so isso? ou o outro intervalo tamb precisa ser considerado?
// }

// double produto(double x, double (*func1)(double), double (*func2)(double))
// {
//     return func1(x) * func2(x);
// }

void decomposicaoLU(double a[MAX], double b[MAX], double c[MAX], double l[MAX], double u[MAX], int n)
{
    u[1] = b[1];
    for (int i = 2; i <= n; i++)
    {
        l[i] = a[i] / u[i - 1];
        u[i] = b[i] - l[i] * c[i - 1];
    }
}

void solucaoLU(double x[MAX], double y[MAX], double d[MAX], double l[MAX], double u[MAX], double c[MAX], int n)
{
    // Ly = d:
    y[1] = d[1];
    for (int i = 2; i <= n; i++)
        y[i] = d[i] - l[i] * y[i - 1];
    // Ux = y;
    x[n] = y[n] / u[n];
    for (int i = n - 1; i > 0; i--)
        x[i] = (y[i] - c[i] * x[i + 1]) / u[i];
}

void imprimir_vetor(double V[MAX], int n, int p)
{
    if (p == 0)
    { // horizontal
        cout << "| ";
        for (int i = 1; i <= n; i++)
        {
            printf("%+7.5f ", V[i]);
            if (i != n)
                cout << " ";
        };
        cout << " |" << endl;
    }
    else
    { // p==1 -> vertical
        for (int i = 1; i <= n; i++)
        {
            cout << " | ";
            printf("%+7.5f ", V[i]);
            cout << "| " << endl;
        }
    }
}