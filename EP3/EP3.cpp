// EP3 MAP3121
// Diogo Vaccaro 8803195
// Otavio Henrique Monteiro 10774159

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>

#define MAX 100
// Constantes
const double PI = 3.141592653589793238463;
const double ks = 3.6; // Condutividade Termica do Silicio W/mK
const double ka = 60;  // Condutividade Termica do Aluminio W/mK
// Variaveis globais de execucao
char questao = '1';   //'t'; // Variavel universal para escolha de questao
bool debug = true;    // Quando true imprime diversos valores intermediarios
bool captura = false; // Quando true e' usada para capturar dados a serem mostrados de cenarios diferentes
// fronteira
double fa = 0.0;
double fb = 0.0;
// dimensoes
double L = 1;      // 20;     // mm
double height = 2; // mm
double P = 30;     // W

using namespace std;

double integra(double a, double b, double xi, double h, double (*func)(double, double, double));
double phi_i(double h, double x, double xi_ant, double xi, double xi_prox);
double k(double x);
double produto_phi_f(double x, double xi, double h);
double produto_phis_normalizados(double x, double xi, double h);
double produto_phis_normalizados_variacao(double x, double xi, double h);
double funcao_escolhida(double x, double fronteira_a, double fronteira_b);
double solucao_exata(double x, double fronteira_a, double fronteira_b);
double Q_gerado_forcante(double x, double L, double sigma);
double Q_gerado(double P, double L, double height);
void solucaoLU(double x[MAX], double y[MAX], double d[MAX], double l[MAX], double u[MAX], double c[MAX], int n);
void decomposicaoLU(double a[MAX], double b[MAX], double c[MAX], double l[MAX], double u[MAX], int n);
void imprimir_vetor(double V[MAX], int n, int p);

int main()
{
    // INICIALIZAR VARIAVEIS
    //  Parametros Gauss (n=2, w=1)
    double T = 1 / sqrt(3); // Abcissa 2 pontos (uma negativa e outra positiva)
    // Matriz Tridiagonal
    int _n[] = {7, 15, 31, 63}; // Dimensao da matriz
    double a[MAX];
    double b[MAX];
    double c[MAX];
    double d[MAX];
    double l[MAX];
    double u[MAX];
    double x[MAX];
    double y[MAX];

    cout.precision(17);

    cout << "Bem-Vindo ao EP3 de MAP3121-2022 \n";
    // cout << "Digite o numero n de linhas e colunas da matriz Anxn: ";
    // cin >> n;

    // EXECUCAO
    for (int iteracao = 0; iteracao < 4; iteracao++)
    {
        int n = _n[iteracao];
        // Montagem da matriz com produtos internos (integrais)
        double h;
        if (questao == 't')
            h = 1 / ((double)n + 1);
        else
            h = L / ((double)n + 1); // Intervalo de [0, L]
                                     // Na validacao q=0 e k=1
                                     // ~f = f +(b-a)*k' - q*(a + (b-a)*x)

        for (int i = 1; i <= n; i++) // montado a cada linha (inclui poucos pontos desnecessarios, com lixo)
        {

            double aux = (double)i;
            double xi = aux * h;
            double xi_ant = (aux - 1) * h;
            double xi_prox = (aux + 1) * h;

            // Para correcao dos intervalos define-se a seguinte consideracao (utilizada localmente)
            // double _xi = xi / L;

            // Construcao da matriz A
                // NOTE: esta montagem considera k(x)=1 e q(x)=0, fazendo com que os produtos internos dos phis sejam simples
                // b[i] = 2 / h;
                // c[i] = -1 / h;
                // a[i + 1] = c[i];
            b[i] = integra(xi_ant / L, xi / L, -1, h, &produto_phis_normalizados) + integra(xi, xi_prox, -1, h, &produto_phis_normalizados); // TODO: conferir se so' o primeiro intervalo de integracao e' normalizado
            c[i] = a[i + 1] = integra(xi / L, xi_prox / L, -1, h, &produto_phis_normalizados_variacao);
            d[i] = integra(xi_ant, xi, xi, h, &produto_phi_f) + integra(xi, xi_prox, xi, h, &produto_phi_f);
        }

        // Opcional, mostrar variaveis construidas
        if (debug && captura)
        {
            imprimir_vetor(a, n + 1, 1);
            cout << "-------" << endl;
            imprimir_vetor(b, n, 1);
            cout << "-------" << endl;
            imprimir_vetor(c, n, 1);
            cout << "-------" << endl;
            imprimir_vetor(d, n, 1);
            cout << "-------" << endl;
        }

        // Decomposicao LU
        decomposicaoLU(a, b, c, l, u, n);

        // Solucao da matriz tridiagonal
        solucaoLU(x, y, d, l, u, c, n);
        // vetor x ? o vetor de alfas
        if (debug)
        {
            cout << "\n---Vetor de Alfas----" << endl;
            imprimir_vetor(x, n, 1);
            cout << "-------" << endl;
        }

        // Resultado exato e comparacoes
        // Validacao
        double erro_max = 0.0; //||u_barra - u||, pegando o maximo
        for (int i = 1; i <= n; i++)
        {
            double v_barra = 0.0;
            double aux_i = (double)i;
            double xi = aux_i * h;
            for (int j = 1; j <= n; j++)
            {
                double aux_j = (double)j;
                double xj = aux_j * h;
                v_barra += x[j] * phi_i(h, xi, xj - h, xj, xj + h);
            }
            // Normalizando
            double u_barra = v_barra + fa + ((fb - fa) * xi / L);
            double u_exato = solucao_exata(xi, fa, fb);

            if (debug)
            {
                cout << "\nPara x=" << xi << " :  ";
                cout << "u" << i << ": encontrado=" << u_barra << "; exato=" << u_exato << "\n";
            }

            erro_max = max(erro_max, u_barra - u_exato);
        }
        cout << "\nErro maximo para n=" << n << " : " << erro_max << endl;

        // Diferentes capturas de valores para testes de Materiais k(x) [alterar na funcao]
        for (int i = 1; i <= n * 10; i++)
        {
            double v_barra = 0.0;
            double aux_i = (double)i;
            double xi = aux_i * h / 10;
            for (int j = 1; j <= n; j++)
            {
                double aux_j = (double)j;
                double xj = aux_j * h;
                v_barra += x[j] * phi_i(h, xi, xj - h, xj, xj + h);
            }
            if (captura)
            {
                cout << v_barra << endl;
            }
        }
    }

    // Finalizar
    cout << "\n\nDigite algum caractere para finalizar.\n";
    char end;
    cin >> end;
    return 1;
}

// Quadratura Gaussiana de 2 pontos
double integra(double a, double b, double xi, double h, double (*func)(double, double, double))
{
    double integral = 0;      // valor da integral
    double ba2 = (b - a) / 2; // valor medio do intervalo ab
    double T = 1 / sqrt(3);   // Abcissa 2 pontos (uma negativa e outra positiva)
    double x, y;
    for (int i = 1; i <= 2; i++) // 2 Pontos
    {
        if (i == 1)
            x = a + ba2 * (-T + 1);
        else
            x = a + ba2 * (T + 1);
        integral += func(x, xi, h);
    }
    integral = integral * ba2; // valor final
    return integral;
}

double phi_i(double h, double x, double xi_ant, double xi, double xi_prox)
{
    double phi;
    if (x <= xi_ant || x > xi_prox)
        phi = 0.0;
    else
    {
        if (x > xi_ant && x <= xi)
            phi = (x - xi_ant) / h;
        else
            phi = (xi_prox - x) / h;
    }
    return phi;
}

double k(double x)
{
    // k=1 (overwriting variacao de material)
    if (questao == '1')
        return 1;

    else if (questao == '2')
        return exp(x);

    // Variacao de material
    double d = L / 4;
    if (x > L / 2 - d && x < L / 2 + d)
    {
        return ks;
    }
    else
    {
        return ka;
    }
}

double produto_phi_f(double x, double xi, double h)
{
    return phi_i(h, x, xi - h, xi, xi + h) * funcao_escolhida(x, fa, fb);
}

double produto_phis_normalizados(double x, double xi, double h)
{
    return k(x) / (h * h * L);
}

double produto_phis_normalizados_variacao(double x, double xi, double h)
{
    return -k(x) / (h * h * L);
}

double Q_gerado(double P, double L, double height) { return P / (L * height); }
double Q_gerado_forcante(double x, double L, double sigma)
{
    double Q0 = 0.5; // W/mm?
    return Q0 * exp(-pow(x - L / 2, 2) / pow(sigma, 2));
}

// FUNCAO PARA CADA QUESTAO
double funcao_escolhida(double x, double fronteira_a, double fronteira_b)
{
    double y;
    switch (questao)
    {
    case '0': // f(x) de Validacao 4.2
        return (12 * x * (1 - x)) - 2;
    case '1': // Para condicoes de fronteira nao homogeneas
        y = (12 * x * (1 - x)) - 2;
        return y + fronteira_a + ((fronteira_b - fronteira_a) * x / L);
    case '2': // exemplo complementar 4.2
        return exp(x)+1;
    case 't': // TESTE
        return x;
    default:
        return 1;
    }
}

double solucao_exata(double x, double fronteira_a, double fronteira_b)
{
    double y;
    switch (questao)
    {
    case '0': // u(x) de Validacao 4.2
        return x * x * (1 - x) * (1 - x);
    case '1': // Para condicoes de fronteira nao homogeneas
        y = x * x * (1 - x) * (1 - x);
    case '2': // exemplo 4.2 complementar
        y = (x-1)*(exp(-x)-1);

        return y + fronteira_a + ((fronteira_b - fronteira_a) * x / L);
    case 't': // TESTE
        return -1;
    default:
        return 1;
    }
}

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
