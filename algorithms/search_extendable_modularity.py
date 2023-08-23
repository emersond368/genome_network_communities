

class SearchExtendableModularity(BaseExtendableHeatmap):


    def greedy_louvain(self,B: np.ndarray, seed: int = None):
        """
        Finds the optimal community structure of the network by maximizing the number of
        within community interactions and minimizing the number of between community interactions.
        This function is a multi-iterative generalization of the louvain community detection
        algorithm. Python implementation is credited to Shi Gui and the group of Dani Bassett.

        Arguments:
        B {np.ndarray} -- NxN Objective function matrix. Typically this is a network adjacency matrix that
        has been corrected by some null nodel. 

        Keyword Arguments:
        seed {[type]} -- A seed that can be used to set the random state for reproducibility purposes (default: {None})

        Raises:
        ValueError: Raised if the maximum iterations to find a community structure is surpassed.

        Returns:
        ci : Nx1 np.array
            final community structure
        q : float
            optimized modularity score
    
        """

        MAX_ATTEMPTS = 1000
        np.random.seed(seed)

        n = len(B)
        ci = np.arange(n) + 1
        Mb = ci.copy()

        B = np.squeeze(np.asarray(B))
        B = (B + B.T)/2.0
        Hnm = np.zeros((n, n))
        for m in range(1, n + 1):
            Hnm[:, m - 1] = np.sum(B[:, ci == m], axis=1)  # node to module degree
        H = np.sum(Hnm, axis=1)  # node degree
        Hm = np.sum(Hnm, axis=0)  # module degree

        q0 = -np.inf
        # compute modularity
        q = np.sum(B[np.tile(ci, (n, 1)) == np.tile(ci, (n, 1)).T])

        first_iteration = True

        while q - q0 > 1e-10:
            it = 0
            flag = True
            while flag:
                it += 1
                if it > MAX_ATTEMPTS:
                    raise ValueError(
                    'Reached maximum attempts to find a partition')
                flag = False
                for u in np.random.permutation(n):
                    ma = Mb[u] - 1
                    dQ = Hnm[u, :] - Hnm[u, ma] + B[u, u]  # algorithm condition
                    dQ[ma] = 0
 
                 max_dq = np.max(dQ)
                 if max_dq > 1e-10:
                     flag = True
                     mb = np.argmax(dQ)
                     Hnm[:, mb] += B[:, u]
                     Hnm[:, ma] -= B[:, u]  # change node-to-module strengths

                     Hm[mb] += H[u]
                     Hm[ma] -= H[u]  # change module strengths

                     Mb[u] = mb + 1

        _, Mb = np.unique(Mb, return_inverse=True)
        Mb += 1

        M0 = ci.copy()
        if first_iteration:
            ci = Mb.copy()
            first_iteration = False
        else:
            for u in range(1, n + 1):
                ci[M0 == u] = Mb[u - 1]  # assign new modules

        n = np.max(Mb)
        b1 = np.zeros((n, n))
        for i in range(1, n + 1):
            for j in range(i, n + 1):
                # pool weights of nodes in same module
                bm = np.sum(B[np.ix_(Mb == i, Mb == j)])
                b1[i - 1, j - 1] = bm
                b1[j - 1, i - 1] = bm
        B = b1.copy()

        Mb = np.arange(1, n + 1)
        Hnm = B.copy()

        q0 = q
        q = np.trace(B)  # compute modularity

    return ci, q

    def leiden(self,B: np.ndarray, seed: int = None):


    return ci, q
