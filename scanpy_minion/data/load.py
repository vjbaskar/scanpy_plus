import anndata
import scanpy as sc

def load_skin(data_name = None):
    data = {
            'fski':'/lustre/scratch126/cellgen/team298/SharedFolders/fetal_skin/data_for_Issac/FS_new_Gene_raw_sliced_processed.h5ad',
            'org':'/lustre/scratch126/cellgen/team298/SharedFolders/fetal_skin/data_for_Issac/Organoid_new/Organoid_new_Gene_raw_sliced_processed_noday13.h5ad',
            'healthy':'/lustre/scratch126/cellgen/team298/SharedFolders/fetal_skin/data_for_Issac/Healthy_all_data.h5ad'
            }
    data_meta = {
            'fski':'fetal skin data',
            'org':'organoid day 13',
            'healthy':'healthy skin'
            }
    if data_name: 
        print(pd.DataFrame(data_meta))

