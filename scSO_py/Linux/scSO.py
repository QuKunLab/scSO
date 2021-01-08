import numpy as np
import pandas as pd
from scipy.sparse.linalg import svds
import warnings
import numpy.matlib
import cvxpy as cp
import scipy.sparse as ss
import libSNMF_inner as SI
import matplotlib.pyplot as plt
import os
eps = 2.22044604925031e-16


class scSO:
    """
    This is a scSo class
    """
    def __init__(self):
        self.Initial_Data = np.array(0)  # save the initial data
        self.Data = np.array(0)
        self.W0 = np.array(0)
        self.n_rank = 2
        self.W = np.array(0)
        self.H = np.array(0)
        self.eigenVector_SM = np.array(0)
        self.eigen_SM = np.array(0)
        self.Fielder_vector = np.array(0)
        self.geneName = np.array(0)
        self.n_spectral = 0

    def read_data(self, data_path):
        """
        This function is used to read csv or tsv file
        """
        if data_path[-3] == 't':
            self.Initial_Data = pd.read_csv(filepath_or_buffer=data_path,
                                            sep='\t',
                                            header=0,
                                            index_col=0)
            geneName = self.Initial_Data.index
            if np.size(geneName) > np.size(np.unique(geneName)):
                geneName = np.unique(geneName)
                self.Initial_Data = self.Initial_Data.groupby(level=0).sum()

            self.Data = self.Initial_Data.values

        elif data_path[-3] == 'c':
            self.Initial_Data = pd.read_csv(filepath_or_buffer=data_path,
                                            sep=',',
                                            header=0,
                                            index_col=0)

            geneName = self.Initial_Data.index
            if np.size(geneName) > np.size(np.unique(geneName)):
                geneName = np.unique(geneName)
                self.Initial_Data = self.Initial_Data.groupby(level=0).sum()

            self.Data = self.Initial_Data.values
        else:
            print("Unsupported file type")

    def filter_genes(self, bound_low=0.1, bound_up=8.5):
        A_ing = self.Initial_Data.values.copy()
        A_ing[A_ing > 0] = 1
        A_ing[A_ing <= 0] = 0
        Row_mean = np.mean(A_ing, axis=1)
        self.Data = self.Initial_Data.values[(Row_mean > 0) &
                                             (Row_mean < 1), :].copy()
        Row_mean = np.mean(self.Data, axis=1)
        rho = np.mean(Row_mean)
        self.Data = self.Data[(Row_mean > bound_low * rho) &
                              (Row_mean < bound_up * rho), :]

    def LogNormalize(self, scale_factor=10000, b_log1p=True):
        """
        normalizes the feature expression measurements
        for each cell by the total expression, multiplies
        this by a scale factor (10,000 by default),
        and log-transforms the result.
        """
        self.Data = np.dot(self.Data,
                           np.diag(scale_factor / np.sum(self.Data, axis=0)))
        self.Data = np.log10(self.Data + 1)

    def calc_W0_and_rank(self, bound=0.085):
        # This Function is used to calculate the initial W and the rank of Data
        W0, single_value, _ = svds(self.Data,
                                   k=min(100,
                                         np.shape(self.Data)[1] - 1),
                                   which='LM')

        single_value = single_value[range(np.shape(W0)[1] - 1, 0, -1)]
        W0 = W0[:, range(np.shape(W0)[1] - 1, 0, -1)]
        W0_P = np.maximum(W0, 0)
        W0_N = np.maximum(-W0, 0)
        n_P_big = np.sum(W0_P * W0_P, axis=0)
        n_N_big = np.sum(W0_N * W0_N, axis=0)
        if np.sum(n_P_big >= n_N_big) > 0:
            W0[:, n_P_big >= n_N_big] = W0_P[:, n_P_big >= n_N_big]

        if np.sum(n_P_big < n_N_big) > 0:
            W0[:, n_P_big < n_N_big] = W0_N[:, n_P_big < n_N_big]

        single_value = single_value[single_value > 0]
        n_svd = min(50, len(single_value) - 1)
        single_value = single_value / single_value[0]
        latent_new_diff = single_value[0:n_svd] / single_value[range(1, 1+n_svd)] - 1
        n_rank = 2
        for i in range(len(latent_new_diff) - 10):
            n_rank = i + 1
            if (latent_new_diff[i] >= bound) and all(latent_new_diff[range(1+i, i+11)] < bound):
                break

        n_rank = max(n_rank, 3)
        self.n_rank = n_rank
        self.W0 = W0

    def lsqnonneg(self, C=None, d=None, tol=None):
        """
        min_{x>0}||Cx-d||_2^2
        Reference: Lawson and Hanson,
        "Solving Least Squares Problems",
        Prentice-Hall, 1974.
        """
        n = np.shape(C)[1]
        # Initialize vector of n zeros and Infs (to be used later)
        nZeros = np.zeros(n)
        wz = nZeros.copy()
        # Initialize set of non-active columns to null
        P = np.zeros(n) == 1
        # Initialize set of active columns to all
        # and the initial point to zeros
        Z = np.ones(n) == 1
        x = nZeros.copy()
        w = d - np.dot(C, x)
        # Set up iteration criterion
        outeriter = 0
        iter = 0
        itmax = 3 * n
        # Outer loop to put variables into set to hold positive coefficients
        while any(Z) and any(w[Z] > tol):
            outeriter = outeriter + 1
            z = nZeros.copy()
            # wz must have the same size as w to preserve
            # the correct indices, so
            # set multipliers to -Inf for variables outside of the zero set.
            wz[P] = float('-Inf')
            wz[Z] = w[Z]
            # Find variable with largest Lagrange multiplier
            t = wz == np.max(wz)
            # Move variable t from zero set to positive set
            P[t] = True
            Z[t] = False
            if any(P):
                if sum(P) < 2:
                    z[P] = d[P] / C[P, P]
                else:
                    C_ing = C[:, P]
                    C_ing = C_ing[P, :]
                    z[P] = np.linalg.solve(C_ing, d[P])

            while any(z[P] <= 0):
                iter = iter + 1
                if iter > itmax:
                    warnings.warn('optimfun:lsqnonneg: IterationCountExceeded')
                    x = z.copy()
                    return x
                # Find indices where intermediate solution z
                # is approximately negative
                Q = (z <= 0) & P
                alpha = np.min(x[Q] / (x[Q] - z[Q]))
                x = x + alpha * (z - x)
                Z = ((np.abs(x) < tol) & P) | Z
                P = ~Z
                z = nZeros.copy()
                if any(P):
                    if sum(P) < 2:
                        z[P] = d[P] / C[P, P]
                    else:
                        C_ing = C[:, P]
                        C_ing = C_ing[P, :]
                        z[P] = np.linalg.solve(C_ing, d[P])

            x = z.copy()
            w = d - np.dot(C, x)
        return x, w

    def SNMF_inner(self, C_mat_=None, B_mat_=None, tol_=None, H=None):
        for i in range(np.shape(H)[1]):
            H[:, i], _ = self.lsqnonneg(C=C_mat_, d=B_mat_[:, i], tol=tol_)

        return H

    def SNMF(self,
             W0=None,
             n_rank=None,
             max_iter=500,
             alpha=0.05,
             beta=0.085,
             bound_error=1e-07):
        '''
        This program is used to implement sparse matrix
        non-negative decomposition algorithm
        '''
        if not W0:
            W0 = self.W0
        if not n_rank:
            n_rank = self.n_rank

        m, n = np.shape(self.Data)
        W = W0[:, 0:n_rank]
        H = np.ones((n_rank, n))
        A_F = np.sum(np.sum(self.Data * self.Data))
        error_new = np.sum(np.sum(H * (np.dot(
            (np.dot(W.T, W)), H)))) - 2 * np.sum(
                np.sum((H * (np.dot(W.T, self.Data)))))
        error_old = 1 + error_new / A_F
        error_min = error_old.copy()
        error_new = error_old.copy()
        gama = 1 - alpha - beta
        iter = 1
        count = 0
        W_old = W.copy()
        H_old = H.copy()
        count_rise = 0
        state = 0
        ones_1_k = np.ones((1, n_rank))
        A_new_H = np.concatenate((gama * self.Data, np.zeros((1, n))), axis=0)
        A_new_W = np.concatenate((gama * self.Data.T, np.zeros((1, m))),
                                 axis=0)
        # print('Iter， error， difference, count, count_rise')
        while count < 5 and iter <= max_iter:

            W_new = np.concatenate((gama * W, beta * ones_1_k), axis=0)
            C_mat = np.dot(W_new.T, W_new)
            B_mat = np.dot(W_new.T, A_new_H)
            tol = 10 * np.max(np.sum(W_new, axis=0)) * (m + 1) * eps
            # self.SNMF_inner(C_mat, B_mat, tol, H)
            H = SI.SNMF_inner(C_mat, B_mat, tol)

            H_new = np.concatenate((gama * H.T, alpha * ones_1_k), axis=0)
            C_mat = np.dot(H_new.T, H_new)
            B_mat = np.dot(H_new.T, A_new_W)
            tol = 10 * np.max(np.sum(H_new)) * n * eps
            # self.SNMF_inner(C_mat, B_mat, tol, W.T).T
            W = SI.SNMF_inner(C_mat, B_mat, tol).T

            error_new = np.sum(H * np.dot(np.dot(W.T, W), H)) - 2 * np.sum(
                H * (np.dot(W.T, self.Data)))

            error_new = 1 + error_new / A_F
            error = abs(error_new - error_old)
            if error_new > (error_old + bound_error):
                count_rise = count_rise + 1
            else:
                count_rise = 0
            if error_new < error_min:
                error_min = error_new.copy()
                W_old = W.copy()
                H_old = H.copy()
            if count_rise > 3:
                print(
                    'arrival the local minimum.'
                )
                error = error_new.copy()
                W = W_old.copy()
                H = H_old.copy()
                return W, H, error_new, state
            if (error <= bound_error):
                count = count + 1
            else:
                count = 0
            # print('Iter %d， %f， %0.5e, %d, %d' %
            #       (iter, error_new, abs(error_new - error_old), count,
            #        count_rise))
            error_old = error_new.copy()

            iter = iter + 1
            if iter == max_iter:
                return W, H, error_new, state

        W = W_old.copy()
        H = H_old.copy()
        state = 1
        return W, H, error_new, state

    def scSO_SNMF(self):
        """
        This function is to call SNMF
        """
        self.calc_W0_and_rank()
        self.W, self.H, _, _ = self.SNMF()

    def zero_norm_inner(self, x_old, A_diff, Lambda):
        # This  function is  used to L0 optimization
        n = np.shape(x_old)[0]
        x = cp.Variable(n)
        y = cp.Variable(n-1)
        cost = (1-Lambda)*cp.norm2(x_old - x) + Lambda*cp.norm1(y)
        prob = cp.Problem(cp.Minimize(cost), [y == A_diff@x])

        prob.solve(solver='SCS')

        return x.value

    def zero_norm_inner_NEW(self, x_old, weight, Lambda):
        # This  function is  used to L0 optimization
        n = np.shape(x_old)[0]
        x = cp.Variable(n)
        y = cp.Variable(n-1)
        cost = (1-Lambda)*cp.norm2(x_old - x) + Lambda*cp.sum(y)
        prob = cp.Problem(cp.Minimize(cost), [y == cp.multiply(x[1:]-x[:-1], weight), y >= 0])

        # cost = cp.norm2(x_old - x)
        # prob = cp.Problem(cp.Minimize(cost), [cp.norm(cp.multiply(x[1:]-x[:-1], weight), 1) <= Lambda])
        prob.solve(solver='SCS')

        return x.value

    def estimate_n_Spector(self, Lambda=0.01):
        """
        This function is used to estimate the number of Spector used in clustering
        """
        n_eigen = min(30, np.shape(self.eigenVector_SM)[1])
        x_ing = self.eigen_SM[:n_eigen]
        x_old_ = x_ing
        n_old = np.shape(x_old_)[0]
        index_i = np.concatenate((range(n_old-1), range(n_old-1)), axis=0)
        index_j = np.concatenate((range(n_old-1), range(1, n_old)), axis=0)
        value_A = np.concatenate((-np.ones(n_old-1), np.ones(n_old-1)), axis=0)
        A_diff = ss.coo_matrix((value_A, (index_i, index_j)))
        A_diff = A_diff.toarray()
        error = 1
        while error > 0.0001:
            y = x_old_[1:] - x_old_[:-1]
            wight = 1/(y+0.00000001)
            wight = np.diag(wight)
            A_diff_wight = np.dot(wight, A_diff)
            x_new = self.zero_norm_inner(x_old_, A_diff_wight, Lambda)
            error = (x_new - x_old_)
            error = np.sqrt(np.sum(error*error))
            x_old_ = x_new

        x_old_ = np.round(x_old_, 3)
        plt.figure(figsize=(6, 6))
        plt.plot(range(np.shape(x_old_)[0]), x_old_, color="blue", linewidth=1.0)
        plt.scatter(range(np.shape(x_old_)[0]), x_ing, s=3, c='k')
        plt.title('Eigen values of Laplace Matrix of cell-cell similarity matrix')
        plt.xlabel('Index of eigen value')
        plt.ylabel('Eigen value')
        plt.show()

        singleValue_bound = x_old_[1]
        # n_eigenVector = np.sum(x_ing < (1.04*singleValue_bound))
        # n_eigenVector = min(n_eigenVector, np.sum(self.eigenVector_SM < 0.01))
        n_eigenVector = np.sum(x_ing < (1.5*singleValue_bound))
        F_v = np.mean(self.eigenVector_SM[:, :n_eigenVector], axis=1)
        if (max(F_v) - min(F_v)) < 1e-4:
            self.eigenVector_SM[:, 0] = 0
            n_eigenVector = n_eigenVector + 1

        n_spectral = min(max(n_eigenVector, 2), self.n_rank-1)
        Fielder_vector = np.mean(self.eigenVector_SM[:, :n_spectral], axis=1)

        x_old = np.sort(Fielder_vector)
        y = x_old[1:] - x_old[:-1]
        block_size = max(int(0.025*np.shape(x_old)[0]), 3)
        for j in range(block_size, np.shape(y)[0] - block_size):
            if y[j] != np.max(y[j-block_size: j + block_size]):
                y[j] = 0
        return n_spectral

    def GMM_fit_and_calcu_0Space(self):
        np.savetxt('H.txt', self.H, delimiter=',')
        main = "./run_run_GMM_MATLAB.sh"
        args_os = " /usr/local/MATLAB/MATLAB_Runtime/v97"
        r_v = os.system(main+args_os)
        if r_v == 0:
            self.eigenVector_SM = pd.read_csv(
                filepath_or_buffer='./EigenVectors.txt',
                sep=',',
                header=None,
                index_col=None).values
            self.eigen_SM = pd.read_csv(
                filepath_or_buffer='./EigenValues.txt',
                sep=',',
                header=None,
                index_col=None).values
            self.eigen_SM = self.eigen_SM[:, 0]
            self.n_spectral = self.estimate_n_Spector(Lambda=0.01)

    def dealWrongCluster(self, x=None, x_old=None):
        x_old = np.round(x_old, 4)
        x = np.round(x, 3)
        y = x.copy()
        x_unique = np.unique(x)
        for i in range(0, np.size(x_unique)):
            if np.size(np.unique(x_old[x == x_unique[i]])) <= 3:
                continue
            x_ing = x_old[x == x_unique[i]]
            length_x = np.sum(x == x_unique[i])
            index_1 = max(int(0.25*length_x), 2)
            index_2 = max(int(0.75*length_x), 4)
            x_ing = x_ing[range(index_1-1, index_2)] - x_ing[index_1-1]
            sig = 2*np.mean(x_ing) / np.size(x_ing)
            if abs(sig) > 0.012:
                if i == 0:
                    y[x == x_unique[i]] = x_unique[i + 1]
                    x_unique[i] = x_unique[i + 1]
                    continue
                if i == np.size(x_unique)-1:
                    y[x == x_unique[i]] = x_unique[i - 1]
                    x_unique[i] = x_unique[i - 1]
                    continue
                if (x_unique[i] - x_unique[i - 1]) > (x_unique[i + 1] - x_unique[i]):
                    y[x == x_unique[i]] = x_unique[i + 1]
                    x_unique[i] = x_unique[i + 1]
                else:
                    y[x == x_unique[i]] = x_unique[i - 1]
                    x_unique[i] = x_unique[i-1]
        return y

    def calcu_FielderVector(self, V_s, Lambda, block_size=None):
        """
        This function is used to calculate  the piecewise constant approximation of V_s
        V_s: is a vector will be approached
        lambda: (scalar)  L0 parameter
        """

        x_old_ = V_s

        if np.shape(V_s)[0]<4000:
            n_old = np.shape(V_s)[0]
            index_i = np.concatenate((range(n_old-1), range(n_old-1)), axis=0)
            index_j = np.concatenate((range(n_old-1), range(1, n_old)), axis=0)
            value_A = np.concatenate((-np.ones(n_old-1), np.ones(n_old-1)), axis=0)
            A_diff = ss.coo_matrix((value_A, (index_i, index_j)))
            A_diff = A_diff.toarray()
        

        error = 1
        while error > 0.0001:
            y = x_old_[1:] - x_old_[:-1]
            if block_size is None:
                block_size = max(np.int(np.ceil(0.01*len(y))), 6)
            y[:block_size] = 0
            y[-block_size:] = 0
            for j in range(block_size-1, np.shape(y)[0]-block_size):
                if y[j] != np.max(y[j-block_size+1:(j+block_size+1)]):
                    y[j] = 0.0

            wight = 1 / (y**1.4+0.0000000001)
            if np.shape(V_s)[0]<4000:
                wight = np.diag(wight)
                A_diff_wight = np.dot(wight, A_diff)
                x_new2 = self.zero_norm_inner(x_old_, A_diff_wight, Lambda)
            else:
                x_new2 = self.zero_norm_inner_NEW(x_old_, wight, Lambda)

            error = np.linalg.norm(x_new2 - x_old_)
            x_old_ = x_new2
        x_new2 = self.dealWrongCluster(x_new2, V_s)
        return x_new2

    def Cluster_func(self, Lambda=0.12):
        """
        This function is used to execute the L0 optimization cluster
        If Lambda>1: cluster based on the top int(Lambda) step clustering
        else L0 optimization clustering
        """
        n_spectral_ing = self.n_spectral

        self.Fielder_vector = np.mean(self.eigenVector_SM[:, :n_spectral_ing], axis=1)
        x_index1 = np.argsort(self.Fielder_vector)
        x_old = self.Fielder_vector[x_index1]
        x_new = x_old.copy()
        state_plus = 0
        if int(Lambda) >= 1:
            if int(Lambda) == 1:
                Lambda = Lambda + 1
            x_old = x_old - np.mean(x_old)
            x_old = x_old / np.abs(np.min(x_old))*2
            y = x_old[1:] - x_old[:-1]
            block_size = max(int(0.025*np.shape(y)[0]), 5)
            for j in range(block_size-1, np.shape(y)[0]-block_size):
                if y[j] != np.max(y[(j-block_size + 1):(j+block_size+1)]):
                    y[j] = 0
            index_y = np.argsort(-y[(block_size-1):-block_size])
            index_X_old = index_y[0:int(Lambda)] + block_size-1
            subBound = 0.5 * (x_old[index_X_old] + x_old[index_X_old+1])
            subBound = np.sort(subBound)
            x_new[x_old < subBound[0]] = np.mean(x_old[x_old < subBound[0]])
            for i in range(1, int(Lambda)-1):
                x_new[x_old >= subBound[i-1] & x_old < subBound[i]] = np.mean(x_old[x_old >= subBound[i-1] & x_old < subBound[i]])

            x_new[x_old >= subBound[int(Lambda) - 1]] = np.mean(x_old[x_old >= subBound[int(Lambda) - 1]])
            x_index = x_index1
        else:
            x_old = x_old - np.mean(x_old)
            x_old = x_old/abs(np.min(x_old))*2
            x_new1 = self.calcu_FielderVector(x_old, Lambda)
            x_new1 = np.round(x_new1, 3)  # Rounding off
            x_new = x_new1
            x_index = x_index1
            # 判断nn_eigenVector+1个最小的奇异值向量的平均作为fielder向量是否可以得到更多的类
            ing = self.eigen_SM[n_spectral_ing] - self.eigen_SM[n_spectral_ing - 1]
            if ((((ing / self.eigen_SM[n_spectral_ing - 1]) < 11)
                 and self.eigen_SM[n_spectral_ing] < 0.02)
                    or (n_spectral_ing < 3)):
            # if (((self.eigen_SM[n_spectral_ing] - self.eigen_SM[n_spectral_ing-1]) < 0.03) and (n_spectral_ing < max(self.n_rank - 1, 3))):
                x_old = np.mean(self.eigenVector_SM[:, :(n_spectral_ing+1)], axis=1)
                x_old = x_old - np.mean(x_old)
                x_old = x_old / abs(min(x_old))*2
                x_index2 = np.argsort(x_old)
                x_old = x_old[x_index2]
                x_new2 = self.calcu_FielderVector(x_old, Lambda)
                x_new2 = np.round(x_new2, 3)
                if np.shape(np.unique(x_new2))[0] > np.shape(np.unique(x_new1))[0]:
                    state_plus = 1
                    x_new = x_new2
                    x_index = x_index2
                    n_spectral_ing = n_spectral_ing + 1

            if state_plus == 0:
                x_old = np.mean(self.eigenVector_SM[:, :(n_spectral_ing-1)], axis=1)
            else:
                x_old = np.mean(self.eigenVector_SM[:, :(n_spectral_ing-2)], axis=1)

            if ((n_spectral_ing - state_plus) > 2 and (np.max(x_old) - np.min(x_old)) > 1e-5):
                x_old = x_old - np.mean(x_old)
                x_old = x_old / np.abs(np.min(x_old))*2
                x_index2 = np.argsort(x_old)
                x_old = x_old[x_index2]
                x_new2 = self.calcu_FielderVector(x_old, Lambda)
                x_new2 = np.round(x_new2, 3)
                if (np.shape(np.unique(x_new2))[0] >= np.shape(np.unique(x_new))[0]):
                    x_new = x_new2
                    x_index = x_index2
                    n_spectral_ing = n_spectral_ing - state_plus - 1

            x_old = np.mean(self.eigenVector_SM[:, :n_spectral_ing], axis=1)
            x_old = x_old - np.mean(x_old)
            x_old = x_old / np.abs(np.min(x_old))*2
            x_index = np.argsort(x_old)

        plt.figure(figsize=(6, 6))
        plt.plot(range(np.shape(x_new)[0]), x_new, color="blue", linewidth=1.0)
        plt.scatter(range(np.shape(x_new)[0]), np.sort(x_old), s=3, c='k')
        plt.title('Eigen Vector of Laplace Matrix of cell-cell similarity matrix')
        plt.xlabel('Index of element in eigen vector')
        plt.ylabel('Value of eigen vector')
        plt.show()

        x_new = np.round(x_new, 3)
        x_new_set = np.unique(x_new)
        celltype_label = np.zeros(np.shape(x_new)[0])
        for i in range(np.shape(x_index)[0]):
            celltype_label[x_index[i]] = np.where(x_new[i] == x_new_set)[0]

        return celltype_label
