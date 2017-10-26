import numpy as np

L = 0.345678
X = 10
T = 10


def compute(Ux, Uxp1, Uxm1, Ut):
    print(2 * (1 - L) * Ux + L * Uxp1 + L * Uxm1 - Ut)


U = np.zeros((X + 1, T + 2), dtype=np.double)
# init
for x in range(1, X):
    U[x][0] = np.sin(x * np.pi / X)
    U[x][1] = np.sin(x * np.pi / X) * np.cos(np.pi / T)

for x in range(1, X):
    for t in range(1, T + 1):
        U[x][t + 1] = 2 * (1 - L) * U[x][t] + L * U[x + 1][t] + L * U[x - 1][t] - U[x][t - 1]

print(U)
print(sum(U[x][T+1] for x in range(X)))
#
# compute(0, 0, 0, 1)
# compute(-1, 0, 0, 0)
#
# compute(0.5, 0, 0, 1)
# compute(-0.952369440632, 0, 0, -0.34567799999999993)
