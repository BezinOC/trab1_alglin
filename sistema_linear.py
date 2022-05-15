import pandas as pd
import numpy as np

class SistemaLinear:


        def __init__(self, n, icod, link_a, link_b, tolm=None):
            self.n = n
            self.icod = icod
            self.a = pd.read_csv(link_a, sep=" ", header=None, dtype=np.float64)
            self.b = pd.read_csv(link_b, sep=" ", header=None, dtype=np.float64)
            self.tolm = tolm


        def foward_subs(self, l, b, n, lu=False):
            y = []
            if lu == True:
                x = 1
            else:
                x = l.loc[0][0]
            y.append(b.loc[0][0] / x)
            for i in range(1, n):
                sum = 0
                for j in range(i):
                    sum += l.loc[i][j] * y[j]
                if lu == True:
                    x = 1
                else:
                    x = l.loc[i][i]
                yi = (b.loc[i][0] - sum) / x
                y.append(yi)
            return pd.DataFrame(y)


        def backward_subs(self, u, y, n):
            x = [0] * n
            x[n-1] = y.loc[n-1][0]/u.loc[n-1][n-1]

            for i in range(n-1,-1,-1):
                sum = 0
                for j in range(i+1, n):
                    sum += u.loc[i][j] * x[j]
                xi = (y.loc[i][0] - sum) / u.loc[i][i]
                x[i] = xi
            return pd.DataFrame(x)


        def decomp_lu(self):
            for k in range(self.n-1):
                for i in range(k+1, self.n):
                    self.a.loc[i][k] = self.a.loc[i][k] / self.a.loc[k][k]
                for j in range(k+1, self.n):
                    for i in range(k+1, self.n):
                        self.a.loc[i][j] = self.a.loc[i][j] - (self.a.loc[i][k]*self.a.loc[k][j])
            y = self.foward_subs(self.a, self.b, self.n, lu=True)
            x = self.backward_subs(self.a, y, self.n)
            return x


        def cholesky(self):
            for i in range(self.n):
                soma = 0
                for k in range(i):
                    soma += (self.a.loc[i][k]) ** 2
                self.a.loc[i][i] = (self.a.loc[i][i] - soma) ** (1/2)
                for j in range(i+1, self.n):
                    soma = 0
                    for k in range(i):
                        soma += self.a.loc[i][k] * self.a.loc[j][k]
                    self.a.loc[j][i] =  (self.a.loc[i][j] - soma) / (self.a.loc[i][i])
                    self.a.loc[i][j] = self.a.loc[j][i]
            y = self.foward_subs(self.a, self.b, self.n)
            x = self.backward_subs(self.a, y, self.n)
            return x


        def jacobi(self):
            x0 = pd.DataFrame([1]*self.n)
            n_interacoes = 0
            err = self.tolm + 1
            while (err > self.tolm):
                x1 = []
                for i in range(self.n):
                    soma = 0
                    for j in range(self.n):
                        if j != i:
                            soma += self.a.loc[i][j] * x0.loc[j][0]
                    x = (self.b.loc[i][0] - soma) / self.a.loc[i][i]
                    x1.append(x)
                x1 = pd.DataFrame(x1)
                norma_x1 = ((x1[0] ** 2).sum() ) ** (1/2)
                norma_x1_x0 = (((x1[0] - x0[0]) ** 2).sum() ) ** (1/2) 
                err = norma_x1_x0 / norma_x1
                n_interacoes += 1
                if n_interacoes > 30:
                    print("Nao convergiu")
                    return -1
                x0 = x1
            print(n_interacoes)
            return x1

        
        def gauss_seidel(self):
            x0 = pd.DataFrame([1.0]*self.n)
            n_interacoes = 0
            err = self.tolm + 1
            while (err > self.tolm):
                norma_x1_x0 = 0
                for i in range(self.n):
                    soma = 0
                    for j in range(self.n):
                        if j != i:
                            soma += self.a.loc[i][j] * x0.loc[j][0]
                    x = (self.b.loc[i][0] - soma) / self.a.loc[i][i]
                    norma_x1_x0 += (x0.loc[i][0] - x)**2
                    x0.loc[i][0] = x
                norma_x0 = ((x0[0] ** 2).sum() ) ** (1/2)
                norma_x1_x0 = norma_x1_x0 ** (1/2) 
                err = norma_x1_x0 / norma_x0
                n_interacoes += 1
                if n_interacoes > 30:
                    print("Nao convergiu")
                    return -1
            print(n_interacoes)
            return x0

# print()
# link_a_lu = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/inputs_sl/A_lu.txt"
# link_b_lu = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/inputs_sl/B_lu.txt"
# s_lu = SistemaLinear(3, 1, link_a_lu, link_b_lu)
# print(s_lu.decomp_lu())
# print()

# link_a_chon = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/inputs_sl/A_cholesky.txt"
# link_b_chon = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/inputs_sl/B_cholesky.txt"
# s_chon = SistemaLinear(3, 1, link_a_chon, link_b_chon)
# print(s_chon.cholesky())
# print()
# #s.test()

# link_a_jacobi = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/inputs_sl/A_jacobi.txt"
# link_b_jacobi = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/inputs_sl/B_jacobi.txt"
# s_jacobi = SistemaLinear(3, 1, link_a_jacobi, link_b_jacobi, 0.001)
# print(s_jacobi.gauss_seidel())