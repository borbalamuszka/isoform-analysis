import numpy as np
import time
from scipy.sparse import csr_matrix, save_npz

n_vars = 60000
nnz = 378257

print("Generating dummy data...")
rows = np.random.randint(0, n_vars, nnz)
cols = np.random.randint(0, n_vars, nnz)
data = np.random.rand(nnz).astype(np.float32)

print("Timing csr_matrix construction...")
t0 = time.time()
mat = csr_matrix((data, (rows, cols)), shape=(n_vars, n_vars), dtype=np.float32)
t1 = time.time()
print(f"csr_matrix took {t1-t0:.4f} seconds.")

print("Timing save_npz...")
t0 = time.time()
save_npz("test_save.npz", mat)
t1 = time.time()
print(f"save_npz took {t1-t0:.4f} seconds.")

import os
os.remove("test_save.npz")
