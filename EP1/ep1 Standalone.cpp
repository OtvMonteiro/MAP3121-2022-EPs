#include<stdio.h>
#include <iostream>
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
    char ciclica;
    int respondido = 0; //bool?
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

    // int n = 4;//20;
    // //lembrar da nomenclatura usada (e.g. d[0] não existe)
    // float a[ ] = { 0, 0, 3, 5, 7};
    // float b[ ] = { 0, 1, 2, 3, 4};
    // float c[ ] = { 0, 1, 3, 5, 0};

    // float d[ ] = { 0, 5, 4, 3, 2};


    // PREENCHER VETORES a,b,c e d
        //Preenchimento pouco eficiente para matrizes grandes
        //Imagino que podemos pegar a entrada de um jeito diferente. Acho que isso facilitaria também, além dos testes, a integração com código dos outros EPs e o código de verificação.
        //Começando com uma estrutura tipo A[MAX][MAX] e rodamos linha a linha (n vezes em cada) para os inputs, tirando depois os paramentros (os termos) das matrizes tridiagonais (ciclicas ou não)
        //Do contrario, imaginando que os parametros serao passados diretamente (as diagonais como vetor) podemos testar com essa estrutura no proprio codigo (sem entrada de usuario): b[ ] = { 0, 1, 3, 5, 7};

        //Acho que podemos pedir no começo se a matriz e' ciclica ou não, e então receber a entrada, alterando-a se necessario. (se possivel, não analisei essa porcao do codigo, mas imagino que ficaria mais elegante)

    /* cout << "Digite o numero n da matriz tridiagonal A nxn: " ;
    cin >> n;

    //Vetor a começa do elemento 1 para manter nomenclatura do pdf
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
    } */

    //Entrando com a matriz inteira
    cout << "Digite o numero n da matriz tridiagonal A nxn: " ;
    cin >> n;
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
    cin >> ciclica;
    if (ciclica == 'S' || ciclica == 's'){ 
        a[1] = A[1][n];
        c[n] = A[n][1];
    }
    else{
        a[1] = 0;
        c[n] = 0;
        //Para calculculo final
        A[1][n] = 0;
        A[n][1] = 0;
    }

    //Entrada de d
    cout << "Digite os valores do termo d: "<<endl; 
    for (int i = 1; i <= n; i++){
        cin >> d[i];
    }

    //Mostra da matriz
    cout <<endl<< "A matriz de entrada e'..."<<endl;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            cout<<A[i][j]<<" ";
        }
        cout << endl;
    }
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
        // cout << "A matriz eh ciclica? [S/N] :"; //anteriormente
        // cin >> ciclica;

        if (ciclica == 'S' || ciclica == 's'){ // resolver matriz c�clica
            float dn = d[n];
            float aux[MAX];
            float v[MAX];
            for(int i=1;i<=n;i++){
                v[i]=0;
            }
            v[1]=a[1];
            v[n]=c[n-1];
            respondido = 1;
            // cout << "Digite o valor do termo a1: ";
            // cin >> a[1];
            // cout << "Digite o valor do termo c" << n << ": ";
            // cin >> c[n];

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

        }else cout <<"Opcao invalida";//Acho desnecessario, da pra considerar só um else acima
    }

    // Imprimir resposta
    cout << "O vetor solucao da matriz A eh:\n";
    imprimir_vetor(x,n); //btw o vetor x é vertical, aqui esta' horizontal

    //Calcular resultado e avaliar se bate
    float solucao[MAX];
    cout<<endl<<endl<<"Comparacao entre solucao encontrada e esperada:"<<endl;
    for (int i = 1; i <= n; i++) {
        solucao[i] = 0;
        for (int j = 1; j <= n; j++) {
            solucao[i] =+ A[i][j]*x[j];
        }
        cout<<"d"<<i<<": encontrada="<<solucao[i]<<"; esperada="<<d[i]<<".\n";
    }


    
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
