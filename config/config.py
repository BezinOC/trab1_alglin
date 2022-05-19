import os

class Settings:

    BASE_PATH = os.getcwd()

    sl_files = {
        1 : {'name': "Decomposição LU", 
            'warning': "AVISO: O método de decomposição LU não funciona para matrizes singulares\n",
            'path': BASE_PATH + "/outputs_sl/out_decomp_lu.txt"},

        2 : {'name': "Cholesky", 
            'warning': "AVISO: O método de Cholesky só funciona para matrizes simétricas positivas definidas\n",
            'path': BASE_PATH + "/outputs_sl/out_cholesky.txt"},

        3 : {'name': "Jacobi", 
            'warning': "AVISO: O método de Jacobi só funciona se A for diagonal dominante\n",
            'path': BASE_PATH + "/outputs_sl/out_jacobi.txt"},

        4 : {'name': "Gauss-Seidel", 
            'warning': "AVISO: O método de Gauss-Seidel só funciona se A for diagonal dominante\n",
            'path': BASE_PATH + "/outputs_sl/out_gauss_seidel.txt"}
    }

    auto_files = {
        1 : {'name': "Power Method", 
            'warning': "AVISO: Ao menos uma componente do vetor inicial deve ser igual a unidade.\n",
            'path': BASE_PATH + "/outputs_auto/power_method.txt"},

        2 : {'name': "Jacobi", 
            'warning': "AVISO: O método de Jacobi só funciona para matrizes simétricas\n",
            'path': BASE_PATH + "/outputs_auto/jacobi.txt"}
    }

    funcao_files = {
        1 : {'name': "Interpolação", 
            'path': BASE_PATH + "/outputs_funcao/interpolacao.txt"},

        2 : {'name': "Regressão", 
            'path': BASE_PATH + "/outputs_funcao/regressao.txt"}
    }

