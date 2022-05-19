import pandas as pd
import numpy as np
from config import config
settings = config.Settings()

class Autovalores:


    def __init__(self, n, link_a, icod, idet, tolm=None):
        self.n = n
        self.a = pd.read_csv(link_a, sep=" ", header=None, dtype=np.float64)
        self.icod = icod
        self.idet = idet
        self.tolm = tolm


    def verificacao(self):
        nlinhas_a, ncolunas_a = len(self.a[0]), len(self.loc.a[0])
        if nlinhas_a != ncolunas_a:
            print("A matriz A não é quadrada") 
            return False
        return True

    
    def write_file(self, autovalores=None, autovetores=None, det=None, n_iteracoes = None):
        with open(settings.auto_files[self.icod]["path"], 'w') as f:
            f.write("Método escolhido: " + settings.auto_files[self.icod]["name"] + "\n")
            f.write(settings.auto_files[self.icod]["warning"])
            if type(autovalores) != pd.DataFrame and type(autovetores) != pd.DataFrame:
                print("Arquivo de saída escrito com erro.")
                return -1
            if self.idet > 0 and det != None:
                f.write("O determinante da matriz é: {}\n".format(det))
            f.write("Número de iterações: {}\n".format(n_iteracoes))
            f.write("Autovalor(es) de A: \n")
            dfAsString = autovalores.to_string(header=False, index=False)
            f.write(dfAsString)
            f.write("\n")
            f.write("Autovetor(es) de A: \n")
            dfAsString = autovetores.to_string(header=False, index=False)
            f.write(dfAsString)
        print("Arquivo de saída escrito com sucesso.")


    def mult_ax(self, a, x):
        y = pd.DataFrame([0.0] * self.n)
        for i in range(self.n):
            y[0] += x.loc[i][0] * self.a[i]
        return y


    def mult_matrix(self, a, b):
        nlinhas_a, ncolunas_a = len(a[0]), len(a.loc[0])
        nlinhas_b, ncolunas_b = len(b[0]), len(b.loc[0])
        if ncolunas_a != nlinhas_b:
            print("Dimensoes erradas")
            return -1
        m = []
        for i in range(nlinhas_a):
            linha = []
            for j in range(ncolunas_b):
                soma = 0
                for k in range(ncolunas_a):
                    soma += a.loc[i][k] * b.loc[k][j]
                linha.append(soma)
            m.append(linha)
        return pd.DataFrame(m)
            

    def power_method(self):
        x0 = pd.DataFrame([1.0] * self.n)
        lamb0 = 1
        n_iteracoes = 0
        err = self.tolm + 1
        while err > self.tolm:
            y = self.mult_matrix(self.a, x0)
            lamb1 = y.max().item()
            x0 =  y / lamb1
            err = abs(lamb1 - lamb0) / abs(lamb1)
            lamb0 = lamb1
            n_iteracoes += 1
            if n_iteracoes > 30:
                self.write_file()
                return -1
        autovalores = pd.DataFrame([lamb0])
        self.write_file(autovalores=autovalores, autovetores=x0, det=None, n_iteracoes=n_iteracoes)
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
                if i != j and abs(m.loc[i][j]) > abs(m.loc[big[0]][big[1]]):
                    big = [i, j]
        return big


    def get_p_matrix(self, big):
        n = self.n
        p = self.identidade(n)
        i, j = big[0], big[1]
        if self.a.loc[i][i] == self.a.loc[j][j]:
            phi = np.pi / 4
        else:
            phi = np.arctan((2*self.a.loc[i][j]) / (self.a.loc[i][i] - self.a.loc[j][j])) / 2
        p.loc[i][i] = np.cos(phi)
        p.loc[j][j] = np.cos(phi)
        p.loc[i][j] = -np.sin(phi)
        p.loc[j][i] = np.sin(phi)
        return p


    def transpose(self, p):
        p_trans = []
        n_columns = len(p.loc[0])
        for i in range(n_columns):
            p_trans.append(p[i])
        return pd.DataFrame(p_trans)


    def jacobi(self):
        x = self.identidade(self.n)
        n_iteracoes = 0
        while(True):
            big = self.get_biggest_out_diag(self.a, self.n)
            if abs(self.a.loc[big[0]][big[1]]) < self.tolm:
                break
            n_iteracoes += 1
            if n_iteracoes > 30:
                self.write_file()
                return -1
            p = self.get_p_matrix(big)
            p_trans = self.transpose(p)
            self.a = self.mult_matrix(p_trans, self.mult_matrix(self.a, p))
            x = self.mult_matrix(x, p)
        autovalores = []
        det = 1
        for i in range(self.n):
            item = self.a.loc[i][i]
            autovalores.append(item)
            det *= item
        autovalores = pd.DataFrame(autovalores)
        if self.idet < 0:
            det = None
        self.write_file(autovalores=autovalores, autovetores=x, det=det, n_iteracoes=n_iteracoes)

    
    def solve(self):
        if self.verificacao == False:
            return -1
        if self.icod == 1:
            return self.power_method()
        elif self.icod == 2:
            return self.jacobi()

        
if __name__ == '__main__':

    filename_a = "A.txt"
    ordem_n = 3
    icod = 2
    idet = 1
    tolm = 0.00001

    link_a = settings.BASE_PATH + "/inputs_auto/" + filename_a

    aut = Autovalores(n = ordem_n, icod = icod, idet = idet, link_a = link_a, tolm=tolm)
    print()
    aut.solve()
    print()


# la = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/a.txt"
# lb = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/b.txt"

# a = pd.read_csv(la, sep=" ", header=None, dtype=np.float64)
# b = pd.read_csv(lb, sep=" ", header=None, dtype=np.float64) 

# print(mult_matrix(a, b))

# link_a = "/Users/bernardocunha/Desktop/2022.1/trab1_alglin/inputs_auto/A.txt"
# a = Autovalores(n=3, link_a=link_a, icod=2, idet=1, tolm = 0.01)
# k = a.jacobi()
#self, n, link_a, icod, idet, tolm=None
