import pandas as pd
import numpy as np
from config import config
settings = config.Settings()

class SistemaLinear:


    def __init__(self, n, icod, idet, link_a, link_b, tolm=None):
        self.n = n
        self.icod = icod
        self.a = pd.read_csv(link_a, sep=" ", header=None, dtype=np.float64)
        self.b = pd.read_csv(link_b, sep=" ", header=None, dtype=np.float64)
        self.tolm = tolm
        self.idet = idet


    def verificacao(self):
        nlinhas_a, ncolunas_a = len(self.a[0]), len(self.loc.a[0])
        nlinhas_b, ncolunas_b = len(self.b[0]), len(self.loc.b[0])
        if nlinhas_a != ncolunas_a:
            print("A matriz A não é quadrada") 
            return False
        elif nlinhas_a != nlinhas_b:
            print("A e B tem número diferente de linhas")
            return False
        elif ncolunas_b != 1:
            print("B não é um vetor coluna")
            return False
        return True


    def write_file(self, x=None, det=None, n_iteracoes = None, hist_err = None):
        with open(settings.sl_files[self.icod]["path"], 'w') as f:
            f.write("Método escolhido: " + settings.sl_files[self.icod]["name"] + "\n")
            f.write(settings.sl_files[self.icod]["warning"])
            if self.idet > 0 and (self.icod == 1 or self.icod == 2):
                f.write("O determinante da matriz é: {}\n".format(det))
            elif self.idet > 0 and (self.icod == 3 or self.icod == 4):
                f.write("O método escohido não permite o cálculo do determinante\n")
            if self.icod == 3 or self.icod == 4:
                f.write("Número de iterações: {}\n".format(n_iteracoes))
                f.write("Histórico da variaçāo do erro: {}\n".format(hist_err))
            if type(x) == pd.DataFrame:
                f.write("O vetor X de saída é: \n")
                dfAsString = x.to_string(header=False, index=False)
                f.write(dfAsString)
                print("Arquivo de saída escrito com sucesso.")
            else:
                print("Arquivo de saída escrito, mas com problemas na execução.")


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


    def is_diag_dom(self, m, n):
        soma_diag = 0
        for i in range(n):
            soma_diag += m.loc[i][i]
        for i in range(n):
            soma_linha, soma_col = 0, 0
            for j in range(n):
                soma_linha += m.loc[i][j]
                soma_col += m.loc[j][i]
            if soma_linha > soma_diag or soma_col > soma_diag:
                return False
        return True


    def decomp_lu(self):
        #decomposição lu
        try:
            for k in range(self.n-1):
                for i in range(k+1, self.n):
                    self.a.loc[i][k] = self.a.loc[i][k] / self.a.loc[k][k]
                for j in range(k+1, self.n):
                    for i in range(k+1, self.n):
                        self.a.loc[i][j] = self.a.loc[i][j] - (self.a.loc[i][k]*self.a.loc[k][j])
            #resolução de x
            y = self.foward_subs(self.a, self.b, self.n, lu=True)
            x = self.backward_subs(self.a, y, self.n)
            #cálculo do det
            if self.idet > 0:
                det = 1
                for i in range(self.n):
                    det *= self.a.loc[i][i]
            else:
                det=None
            #escrita do arquivo
            self.write_file(x=x, det=det)
        except:
            self.write_file()


    def cholesky(self):
        #decomposição
        try:
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
            #resolução de x
            y = self.foward_subs(self.a, self.b, self.n)
            x = self.backward_subs(self.a, y, self.n)
            #cálculo do det
            if self.idet > 0:
                det = 1
                for i in range(self.n):
                    det *= self.a.loc[i][i]
                det = det**2
            else:
                det = None
            #escrita do arquivo
            self.write_file(x=x, det=det)
        except:
            self.write_file()


    def jacobi(self):
        if self.is_diag_dom(self.a, self.n) == False:
            self.write_file()
            return -1
        x0 = pd.DataFrame([1]*self.n)
        n_iteracoes = 0
        err = self.tolm + 1
        hist_err = []
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
            hist_err.append(err)
            n_iteracoes += 1
            if n_iteracoes > 30:
                self.write_file()
                return -2
            x0 = x1
        self.write_file(x1, n_iteracoes=n_iteracoes, hist_err=hist_err)

        
    def gauss_seidel(self):
        if self.is_diag_dom(self.a, self.n) == False:
            self.write_file()
            return -1
        x0 = pd.DataFrame([1.0]*self.n)
        n_iteracoes = 0
        err = self.tolm + 1
        hist_err = []
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
            hist_err.append(err)
            n_iteracoes += 1
            if n_iteracoes > 30:
                self.write_file()
                return -2
        self.write_file(x0, n_iteracoes=n_iteracoes, hist_err=hist_err)

        
    def solve(self):
        if self.verificacao == False:
            return -1
        if self.icod == 1:
            return self.decomp_lu()
        elif self.icod == 2:
            return self.cholesky()
        elif self.icod == 3:
            return self.jacobi()
        elif self.icod == 4:
            return self.gauss_seidel()


if __name__ == '__main__':

    filename_a = "A_cholesky.txt"
    filename_b = "B_cholesky.txt"
    ordem_n = 3
    icod = 2
    idet = 1
    tolm = None

    link_a = settings.BASE_PATH + "/inputs_sl/" + filename_a
    link_b = settings.BASE_PATH + "/inputs_sl/" + filename_b

    sl = SistemaLinear(n = ordem_n, icod = icod, idet = idet, link_a = link_a, link_b = link_b, tolm=tolm)
    print()
    sl.solve()
    print()