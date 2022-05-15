import pandas as pd
import numpy as np

class Autovalores:


        def __init__(self, n, link_a, icod, idet, tolm=None):
            self.n = n
            self.a = pd.read_csv(link_a, sep=" ", header=None, dtype=np.float64)
            self.icod = icod
            self.idet = idet
            self.tolm = tolm

        def mult_ax(self, a, x):
            y = pd.DataFrame([0.0] * self.n)
            for i in range(self.n):
                y[0] += x.loc[i][0] * self.a[i]
            return y


        def power_method(self):
            x0 = pd.DataFrame([1.0] * self.n)
            lamb0 = 1
            n_interacoes = 0
            err = self.tolm + 1
            while err > self.tolm:
                y = self.mult_ax(self.a, x0)
                lamb1 = y.max().item()
                x0 =  y / lamb1
                err = abs(lamb1 - lamb0) / abs(lamb1)
                lamb0 = lamb1
                n_interacoes += 1
                if n_interacoes > 30:
                    print("Falhou")
                    return -1
            print(n_interacoes)
            return x0

        
        def identidade(self, n):
            m = []
            for i in range(n):
                m_linha = []
                for j in range(n):
                    if i == j: 
                        m_linha.append(1.0)
                    else:
                        m_linha.append(0.0)
                m.append(m_linha)
            m = pd.DataFrame(m)
            return m


        def get_biggest_out_diag(self, m, n):
            #m = self.a
            big = [0, 1]
            for i in range(n):
                for j in range(n):
                    if i != j and m.loc[i][j] > m.loc[big[0]][big[1]]:
                        big = [i, j]
            return big


        def get_p_matrix(self, big):
            n = self.n
            p = self.identidade(n)
            i, j = big[0], big[1]
            if self.a.loc[i][j] == self.a.loc[j][i]:
                phi = np.pi / 4
            else:
                phi = np.arctan((2*self.a.loc[i][j]) / (self.a.loc[i][i], self.a.loc[j][j])) / 2
            p.loc[i][i] = np.cos(phi)
            p.loc[j][j] = np.cos(phi)
            p.loc[i][j] = -np.sin(phi)
            p.loc[j][i] = np.sin(phi)
            return p


        def jacobi(self):
            x = self.identidade(self.n)
            k = 1
            while(True):
                big = self.get_biggest_out_diag(self.a, self.n)
                return self.get_p_matrix(big)


            
            
link_a = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/inputs_auto/A.txt"
a = Autovalores(3, link_a, 1, 1, 0.00001)
print(a.jacobi())