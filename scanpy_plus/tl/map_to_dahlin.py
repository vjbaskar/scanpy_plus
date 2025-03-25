from ..globimport import *
import anndata
import sklearn
import matplotlib.pyplot as plt
import sys    


class Project_landscape:
    """
    Projecting onto landscapes and performing label transfer
    For Niki's dataset scale together as they are both 10x
    For Nestorowa's dataset scale separately

    Run methods flowchart:
    -----------------------
    normalise, log, scale_(together/separately), 
    run_pca_ref, pca_project_myadata, calc_pcadist,
    labelling, plot_on_ref, cluster_assignment
    """
    def __init__(self, my_adata, ref_adata, genes_to_consider, ref_data_obs_label = 'CellSubType', npcs = 25, knn = 15):
        OLG = list(set(np.intersect1d(genes_to_consider, my_adata.var_names)))
        print("Genes common to both:")
        print(len(OLG))
        print(OLG[0:4])
        my_adata = my_adata[:,OLG].copy()
        ref_adata = ref_adata[:,OLG].copy()
        self.my_adata = my_adata
        self.ref_adata = ref_adata
        self.genes_to_consider = genes_to_consider
        self.npcs = npcs
        self.dist_mtx = dict()
        self.ref_data_obs_label = ref_data_obs_label
        self.knn = knn

        plt.rcParams["figure.figsize"] = (7,4)
        self.fig = plt.figure()


    def normalise(self):
        """
        Normalise per cell
        """
        sc.pp.normalize_total(self.my_adata, target_sum=1e4, exclude_highly_expressed=True)
        sc.pp.normalize_total(self.ref_adata, target_sum=1e4, exclude_highly_expressed=True)

    def log(self):
        """
        Logp data
        """

        sc.pp.log1p(self.my_adata)
        sc.pp.log1p(self.ref_adata)



    def scale_together(self):
        """
        Combine the two anndata.
        Scale them together. 
        Adviced for Niki's dataset integration
        """
        # scale them together
        data_comb = self.my_adata.concatenate(self.ref_adata)
        sc.pp.scale(data_comb)
        # Split them back again
        self.my_adata = anndata.AnnData(X=data_comb[data_comb.obs['batch'] == '0',:].X, 
                            obs=self.my_adata.obs, 
                            var=self.my_adata.var)
        self.ref_adata = anndata.AnnData(X=data_comb[data_comb.obs['batch'] == '1',:].X, 
                                 obs=self.ref_adata.obs, 
                                 var=self.ref_adata.var, 
                                 obsm=self.ref_adata.obsm, 
                                 uns=self.ref_adata.uns)
        self.data_comb = data_comb

    def scale_separately(self):
        """
        Scales the two data separately.
        Adviced if the datasets are from two different sequencing methods, for eg Nesterowa and 10x
        """
        # Scale them separately because they are from different technologies
        sc.pp.scale(self.my_adata)
        sc.pp.scale(self.ref_adata)


    def run_pca_ref(self):
        """
        Run PCA and optain pca object for ref data
        """
        # PCs for Niki's data
        pca_ = sklearn.decomposition.PCA(n_components=50, svd_solver='auto', random_state=0)
        pca_.fit(self.ref_adata.X)
        # Plot variance contrib
        ax1 = self.fig.add_subplot(121)
        ax1.plot(pca_.explained_variance_)
        self.pca = pca_

    def pca_project_myadata(self):
        """
        Using PCs from ref data predict the coords for my_adata
        """
        # Project data and niki's data onto Niki's PCs
        pca_ = self.pca
        self.my_adata_proj = pca_.transform(self.my_adata.X)
        self.ref_adata_proj = pca_.transform(self.ref_adata.X)
        # Plot projection
        #ax2 = self.fig.add_subplot(122)
        #ax2.scatter(self.ref_adata_proj[:,0], self.ref_adata_proj[:,1], c='black', alpha=0.5)
        #plt.xlabel('PCA1')
        #plt.ylabel('PCA2')
        #ax2.scatter(self.my_adata_proj[:,0], self.my_adata_proj[:,1], c='red', alpha=0.5)
        #plt.show()
        plt.scatter(self.ref_adata_proj[:,0], self.ref_adata_proj[:,1], c='black', alpha=0.5)
        plt.xlabel('PCA1')
        plt.ylabel('PCA2')
        plt.scatter(self.my_adata_proj[:,0], self.my_adata_proj[:,1], c='red', alpha=0.5)

    def calc_pca_dist(self):
        """
        Compute euclidean dist between each cell in my_adata to each cell in ref_adata in the PCA space
        """
        # High memory process
        npcs=self.npcs
        from sklearn.metrics.pairwise import euclidean_distances
        self.pca_dist = euclidean_distances(self.my_adata_proj[:,0:npcs], self.ref_adata_proj[:,0:npcs])

    def labelling(self):
        """
        Using the distances find knn nearest neighbours of ref_adata for each cell in my_adata.
        The most common NN in ref_adata is returned as the label for each cell in my_adata.
        """
        def progressbar(it, prefix="", size=60, out=sys.stdout): # Python3.3+
            count = len(it)
            def show(j):
                x = int(size*j/count)
                print("{}[{}{}] {}/{}".format(prefix, "#"*x, "."*(size-x), j, count), 
                        end='\r', file=out, flush=True)
            show(0)
            for i, item in enumerate(it):
                yield item
                show(i+1)
            print("\n", flush=True, file=out)
        
        knn = self.knn
        from collections import Counter
        from collections import defaultdict
        cl_assigned = []
        D_sub = self.pca_dist
        ref_label = self.ref_data_obs_label
        Rstore = defaultdict(list) # dictionary to store results

        for i in progressbar(range(D_sub.shape[0]), prefix = "Label transfer "):
            CellDis = D_sub[i,:]
            CellDis_sorted = np.argsort(CellDis)[:knn]
            max_samples = self.ref_adata.obs_names[CellDis_sorted]
            cl_assigned.append(max_samples)
            Rstore['MinDist'].append(np.min(CellDis[CellDis_sorted]))
            Rstore['MedianDist'].append(np.median(CellDis[CellDis_sorted]))
            Rstore['MaxDist'].append(np.max(CellDis[CellDis_sorted]))
            Rstore['SD'].append(np.std(CellDis[CellDis_sorted]))
            Rstore['label_ct'].append(Counter(self.ref_adata.obs[ref_label][CellDis_sorted]).most_common(1)[0][0])
        Rstore = pd.DataFrame.from_dict(Rstore)
        Rstore.index = self.my_adata.obs_names
        self.label_store = Rstore
        self.cl_assigned = cl_assigned

    def plot_on_ref(self, my_adata_cluster_obs = 'leiden', prefix = 'niki'):
        """
        Using the distances map each cluster in my_adata onto cells in niki dataset.
        """
        from collections import Counter
        from collections import defaultdict
        
        ref_data = self.ref_adata
        proj_data = self.my_adata
        proj_data_obs = my_adata_cluster_obs
        cl_assigned = self.cl_assigned
        self.my_adata_cluster_obs = my_adata_cluster_obs
        obs_names = []
        #CT = np.unique(proj_data.obs[proj_data_obs])
        #CT = np.sort(proj_data.obs[proj_data_obs].unique().astype('int')) 
        CT = proj_data.obs[proj_data_obs].unique()
        #print(CT)
        for ct in CT:
            ct = str(ct)
            cl_assigned_sub = [cl_assigned[i] for i in np.where(proj_data.obs[proj_data_obs] == ct)[0]]
            cl_flat_sub = [item for sublist in cl_assigned_sub for item in sublist]
            freq1 = Counter(cl_flat_sub)
            freq2 = np.array([0] * ref_data.X.shape[0])
            for k, i in freq1.items():
                if k in ref_data.obs_names:
                    idx = np.where(ref_data.obs_names==k)
                    freq2[idx] = i

            ref_data.obs[prefix+'_'+ct] = np.log2(freq2+1)
            obs_names.append(prefix+'_'+ct)
        self.ref_obsnames = obs_names
        sc.pl.draw_graph(self.ref_adata, color=obs_names, size = 30)

    def cluster_assignment(self):
        """
        Based on cell labelling, use the most common label for the cells in each cluster of my_adata as the cluster name.
        """
        cluster_assignment = dict()
        cluster_obs = self.my_adata_cluster_obs
        Rstore = self.label_store
        ref_adata = self.ref_adata
        my_adata = self.my_adata 
        for i in my_adata.obs[cluster_obs].unique():
            #print(Counter(np.array(cl_assigned_ct)[adata.obs['louvain_v2'] == i]))
            ct = Counter(np.array(Rstore['label_ct'])[my_adata.obs[cluster_obs] == i]).most_common(1)[0][0]
            print(i+':'+ct)
            cluster_assignment[int(i)] = ct
        self.classign = cluster_assignment
 
    
def map_to_dahlin(data):
    
    """
    Maps to dahlin landscape
    Expects the anndata to have anndata.raw instantiated. 
    
    Run methods flowchart:
    -----------------------
        normalise, log, scale_(together/separately), 
        run_pca_ref, pca_project_myadata, calc_pcadist,
        labelling, plot_on_ref, cluster_assignment
    
    Returns: Project_landscape object
    
    """
    
    import sys
    def progressbar(it, prefix="", size=60, out=sys.stdout): # Python3.3+
        count = len(it)
        def show(j):
            x = int(size*j/count)
            print("{}[{}{}] {}/{}".format(prefix, "#"*x, "."*(size-x), j, count), 
                    end='\r', file=out, flush=True)
        show(0)
        for i, item in enumerate(it):
            yield item
            show(i+1)
        print("\n", flush=True, file=out)
    

    # Query landscape
    #adata_query = anndata.AnnData(X=np.exp(data.raw.X.todense())-1,
    #                        obs=data.obs, 
    #                        var=data.raw.var, 
    #                        obsm=data.obsm,
    #                        uns = data.uns)
    adata_query = data.copy()
    adata_query = adata_query.raw.to_adata()

    # niki's data
    adata_ref = sc.read("/rds/project/rds-SDzz0CATGms/users/vs401/refdata/10x/mm/niki_landscape/ref_landscape.h5ad")
    adata_ref.var_names_make_unique()
    niki_hvg = np.genfromtxt('/rds/project/rds-SDzz0CATGms/users/vs401/refdata/10x/mm/niki_landscape/gene_names.txt',dtype=str)
    pob = Project_landscape(my_adata = adata_query.copy(), ref_adata = adata_ref.copy(), genes_to_consider = adata_ref.uns['hvg'])
    pob.normalise()
    pob.log()
    pob.scale_together()
    pob.run_pca_ref()
    pob.pca_project_myadata()
    pob.calc_pca_dist()
    pob.labelling()
    temp = pob.label_store.reindex(data.obs.index)
    data.obs['celltype'] = temp['label_ct']
    return pob
