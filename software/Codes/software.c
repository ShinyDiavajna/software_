 #include <stdio.h>
#include <math.h>

void mul(double X[10][10], double Y[10][10], double Z[10][10], int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Z[i][j] = 0;
            for (int k = 0; k < N; k++) {
                Z[i][j] += X[i][k] * Y[k][j];
            }
        }
    }
}

void gramSchmidt(double X[10][10], double Q[10][10], double R[10][10], int N) {
    double temp[10];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Q[i][j] = 0;
            R[i][j] = 0;
        }
    }

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            temp[i] = X[i][j];
        }

        for (int i = 0; i < j; i++) {
            R[i][j] = 0;
            for (int k = 0; k < N; k++) {
                R[i][j] += Q[k][i] * X[k][j];
            }

            for (int k = 0; k < N; k++) {
                temp[k] -= R[i][j] * Q[k][i];
            }
        }

        R[j][j] = 0;
        for (int k = 0; k < N; k++) {
            R[j][j] += temp[k] * temp[k];
        }
        R[j][j] = sqrt(R[j][j]);

        for (int k = 0; k < N; k++) {
            Q[k][j] = temp[k] / R[j][j];
        }
    }
}

void qrAlgorithm(double X[10][10], double eigenvalues[10], int N) {
    double Q[10][10], R[10][10], temp[10][10];
    double tolerance = 1e-10;
    double diff;

    for (int i = 0; i < N; i++) {
        eigenvalues[i] = X[i][i];
    }

    for (int iteration = 0; iteration < 1000; iteration++) {
        gramSchmidt(X, Q, R, N);

        mul(R, Q, temp, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                X[i][j] = temp[i][j];
            }
        }

        diff = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    diff += X[i][j] * X[i][j];
                }
            }
        }

        if (diff < tolerance) {
            break;
        }
    }

    for (int i = 0; i < N; i++) {
        eigenvalues[i] = X[i][i];
    }
}

int main() {
    int N;

    printf("Enter the size of the matrix (N x N): ");
    scanf("%d", &N);

    double X[10][10], eigenvalues[10];

    
    printf("Enter the elements of the matrix (row by row):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("Enter element X[%d][%d]: ", i + 1, j + 1);
            scanf("%lf", &X[i][j]);
        }
    }

    qrAlgorithm(X, eigenvalues, N);

    printf("Eigenvalues of the matrix are:\n");
    for (int i = 0; i < N; i++) {
        printf("%lf\n", eigenvalues[i]);
    }

    return 0;
}

