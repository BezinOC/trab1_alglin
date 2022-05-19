import pandas as pd
import numpy as np
from config import config
settings = config.Settings()

class FuncaoEmUmPonto:


        def __init__(self, n, icod, link_a, link_b, interpol_val):
            self.n = n
            self.icod = icod
            self.a = pd.read_csv(link_a, sep=" ", header=None, dtype=np.float64)
            self.b = pd.read_csv(link_b, sep=" ", header=None, dtype=np.float64)
            self.interpol_val = interpol_val


        def write_file(self, x, y):
            with open(settings.funcao_files[self.icod]["path"], 'w') as f:
                f.write("Método escolhido: " + settings.funcao_files[self.icod]["name"] + "\n")
                f.write("Para x = {}, temos P({}) = {}".format(x, x, y))
                print("Arquivo de saída escrito com sucesso.")

        def lagrange(self):
            coef = []
            for i in range(self.n):
                L = 1
                for j in range(len(self.a)):
                    if i != j:
                        L = L * (self.interpol_val - self.a.loc[j]) / (self.a.loc[i] - self.a.loc[j])
                coef.append(L)
            
            pn = 0
            for i in range(self.n):
                pn = pn + self.b.loc[i] * coef[i]
            
            #return f'p({self.interpol_val}) = {pn[0]}'
            self.write_file(x = self.interpol_val, y = pn[0])

        def regression(self):
            soma_A = 0
            for i in range(self.n):
                soma_A = soma_A + self.a.loc[i] 
            media_A = soma_A / self.n
            soma_B = 0
            for i in range(self.n):
                soma_B = soma_B + self.b.loc[i] 
            media_B = soma_B / self.n
            error_A = self.a - media_A
            error_B = self.b - media_B
            error_sum = np.sum(error_A * error_B)
            error_A_2 = error_A ** 2
            error_A_2_sum = np.sum(error_A_2)
            m = error_sum / error_A_2_sum
            c = media_B - m * media_A

            result = c[0] + self.interpol_val * m[0]
            #return result 
            self.write_file(x = self.interpol_val, y = result)
        
        def solve(self):
            if self.icod == 1:
                return self.lagrange()
            elif self.icod == 2:
                return self.regression()


if __name__ == '__main__':

    filename_a = "X.txt"
    filename_b = "Y.txt"
    n_pares = 3
    icod = 1
    interpol_val = 2.1

    link_a = settings.BASE_PATH + "/inputs_funcao/" + filename_a
    link_b = settings.BASE_PATH + "/inputs_funcao/" + filename_b

    fun = FuncaoEmUmPonto(n = n_pares, icod = icod, link_a = link_a, link_b = link_b, interpol_val=interpol_val)
    print()
    fun.solve()
    print()

# filename_a = "X.txt"
# filename_b = "Y.txt"

# link_a = settings.BASE_PATH + '/inputs_funcao/' + filename_a
# link_b = settings.BASE_PATH + '/inputs_funcao/' + filename_b 
# func = FuncaoEmUmPonto(3, 2, link_a, link_b, 2.1)
# result_lagr = func.lagrange()
# print(result_lagr)
# print()
# result_regr = func.regression()
# print(result_regr)
