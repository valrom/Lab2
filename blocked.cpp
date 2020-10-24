#include <iostream>
#include <cstdlib>
#include <mpi.h>

using namespace std;

size_t N;

void add (float* arr1, float* arr2, float* rslt);
void sub (float* arr1, float* arr2, float* rslt);
void mul (float* arr1, float* arr2, float* rslt);
void div (float* arr1, float* arr2, float* rslt);

float *A, *B, *C, *D, *E, *G;
float *T0, *T1, *T2;
float *Y1, *Y2, *Y3;

int main(int argc, char **argv) {

    int rank;
    double time0, time1;
    MPI_Status status;

    if ( argc < 2 ) {
        cout << "No argumets" << endl;
        return 0;
    }
    // Количество данных в массивах
    N = atoi( argv[ 1 ] );


    A = new float[N];
    B = new float[N];
    C = new float[N];
    D = new float[N];
    E = new float[N];
    G = new float[N];

    T0 = new float[N];
    T1 = new float[N];
    T2 = new float[N];

    Y1 = new float[N];
    Y2 = new float[N];
    Y3 = new float[N];

    // Инциализация массивов
    for (int i=0;i<N;i++) {
        T0[i]=0;
        T1[i]=0;
        T2[i]=0;

        Y1[i]=0;
        Y2[i]=0;
        Y3[i]=0;

        A[i]=i;
        B[i]=i+1;
        C[i]=2*i+1;
        D[i]=3*i+1;
        E[i]=2*i+2;
        G[i]=3*i+2;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        time0 = MPI_Wtime();
        sub(A, E, T0);//T0 = A-E
        MPI_Send(T0, N, MPI_FLOAT, 1, 1, MPI_COMM_WORLD);//отправка A-E процессу 1
        MPI_Send(T0, N, MPI_FLOAT, 2, 1, MPI_COMM_WORLD);//отправка A-E процессу 2
        MPI_Recv(T1, N, MPI_FLOAT, 1, 1, MPI_COMM_WORLD, &status);//получить в T1 C*G от процесса 1
        add(T0,T1,T0);//T0 = A-E + C*G
        mul(T0,T0,T0);//T0 = (A-E + C*G)^2
        MPI_Recv(T1, N, MPI_FLOAT, 2, 1, MPI_COMM_WORLD, &status);//получить в T1 G^2 от процесса 2
        div(T0, T1, Y1);//Y1 = (A-E + C*G)^2 / G^2
        MPI_Send(Y1, N, MPI_FLOAT, 1, 1, MPI_COMM_WORLD);//отправка Y1 процессу 1
        MPI_Recv(Y2, N, MPI_FLOAT, 1, 1, MPI_COMM_WORLD, &status);//получить в Y2 Y2 от процесса 1
        MPI_Recv(Y3, N, MPI_FLOAT, 2, 1, MPI_COMM_WORLD, &status);//получить в Y3 Y3 от процесса 2
        time1 = MPI_Wtime();
        cout<<"Thread "<<rank<<" with runtime "<<time1-time0<<endl;
    }

    else if(rank == 1){
        time0 = MPI_Wtime();
        mul(C, G, T0);//T0 = C*G
        MPI_Recv(T1, N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);//получить в T1 A-E от процесса 0
        MPI_Send(T0, N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);//отправка C*G процессу 0
        MPI_Send(T0, N, MPI_FLOAT, 2, 1, MPI_COMM_WORLD);//отправка C*G процессу 2
        mul(T1, G, T1);//T1 = (A-E)*G
        sub(T1, C, T1);//T1 = (A-E)*G - C
        MPI_Recv(T0, N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);//получить в T0 Y1 от процесса 0
        sub(T0, T1, Y2);//Y2 = Y1 - [(A-E)*G - C]
        MPI_Send(Y2, N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);//отправка Y2 процессу 0
        time1 = MPI_Wtime();
        cout<<"Thread "<<rank<<" with runtime "<<time1-time0<<endl;
    }

    else if(rank == 2){
        time0 = MPI_Wtime();
        mul(G, G, T0);//T0 = G*G
        MPI_Recv(T1, N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);//получить в T1 A-E от процесса 0
        MPI_Send(T0, N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);//отправка G*G процессу 0
        mul(T1, T0, T0);//T0 = (A-E)*G^2
        MPI_Recv(T2, N, MPI_FLOAT, 1, 1, MPI_COMM_WORLD, &status);//получить в T2 C*G от процесса 1
        add(T2,T0,T0);//T0 = C*G+(A-E)*G^2
        mul(T1,T0,T0);//T0 = [A-E]*[C*G+(A-E)*G^2]
        MPI_Recv(T1, N, MPI_FLOAT, 3, 1, MPI_COMM_WORLD, &status);//получить в T1 B*D*[B*D-E-G]*[B*D-E]^2 от процесса 3
        mul(T0, T1, Y3);//Y3 = [A-E]*[C*G+(A-E)*G^2]*B*D*[B*D-E-G]*[B*D-E]^2
        MPI_Send(Y3, N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);//отправка Y3 процессу 0
        time1 = MPI_Wtime();
        cout<<"Thread "<<rank<<" with runtime "<<time1-time0<<endl;
    }

    else if(rank == 3){
        time0 = MPI_Wtime();
        mul(B, D, T0);//T0 = B*D
        sub(T0, E, T1);//T1 = B*D - E
        sub(T1, G, T2);//T2 = B*D - E - G
        mul(T1, T1, T1);//T1 = [B*D-E]^2
        mul(T1, T2, T1);//T1 = [B*D-E-G]*[B*D-E]^2
        mul(T1, T0, T1);//T1 = B*D*[B*D-E-G]*[B*D-E]^2
        MPI_Send(T1, N, MPI_FLOAT, 2, 1, MPI_COMM_WORLD);//отправка B*D*[B*D-E-G]*[B*D-E]^2  процессу 2
        time1 = MPI_Wtime();
        cout<<"Thread "<<rank<<" with runtime "<<time1-time0<<endl;
    }

    MPI_Finalize();

    // Освобождаем всю выделенную память
    delete []  A;
    delete []  B;
    delete []  C;
    delete []  D;
    delete []  E;
    delete []  G;
    delete [] T0;
    delete [] T1;
    delete [] T2;
    delete [] Y1;
    delete [] Y2;
    delete [] Y3;


    return 0;
}

/*Функции для выполнения операций над матрицами*/

void add (float* arr1, float* arr2, float* rslt)
{
    for (int i=0;i<N;i++)
        rslt[i]= arr1[i]+arr2[i];
}

void sub (float* arr1, float* arr2, float* rslt)
{
    for (int i=0;i<N;i++)
        rslt[i]= arr1[i]-arr2[i];
}

void mul (float* arr1, float* arr2, float* rslt)
{
    for (int i=0;i<N;i++)
        rslt[i]= arr1[i]*arr2[i];
}

void div (float* arr1, float* arr2, float* rslt)
{
    for (int i=0;i<N;i++)
        if (arr2[i]) rslt[i]= arr1[i]/arr2[i];
}