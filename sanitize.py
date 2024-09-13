# import numpy as np
# import pandas as pd
# edgelist = pd.read_csv(PATH_HERE, sep = "\t", header=None)
# edgelist = np.array(edgelist)
# inverse = np.zeros(int(edgelist[0,0]))
# mappings = np.ones(int(edgelist[1:].max()) + 1)*-1
# sanitized = []
# sanitized.append([edgelist[0,0], edgelist[0,1]])
# counter = 0
# for i in tqdm(range(1, len(edgelist))):
#     e = edgelist[i]
#     u = int(e[0])
#     v = int(e[1])
#     if mappings[u] == -1:
#         mappings[u] = counter
#         inverse[counter] = u
#         counter += 1
#     if mappings[v] == -1:
#         mappings[v] = counter
#         inverse[counter] = v
#         counter += 1
#     x = int(mappings[u])
#     y = int(mappings[v])
#     sanitized.append([x, y])

# coo = np.array(sanitized)
# coo_u =  np.sort(coo[1:], axis = 1)
# coo_u = np.unique(coo_u, axis = 0)
# row_indices = np.where(coo_u[:, 0] == coo_u[:, 1])[0]
# if (len(row_indices)):
#     coo_u = np.delete(coo_u, row_indices, axis = 0)
# edges = len(coo_u)
# print(edges, np.max(coo_u)+1)
# vertices = np.max(coo_u) + 1
# coo = np.vstack(([vertices, edges], coo_u))
# np.save(SAVE_MAPPING_TO_OG_COO, inverse)
# np.save(SAVE_RELABELLED_COO, coo)


import numpy as np
import pandas as pd
from tqdm import tqdm

# Set your file paths
INPUT_PATH = "SPANNER/snap-soc-livejournal.txt"  # Path to your input file
OUTPUT_PATH = "SPANNER/snap-soc-livejournal-sanitized.txt"  # Path to save the sanitized file
SAVE_MAPPING_TO_OG_COO = "SPANNER/mapping_to_og.npy"  # Path to save the mapping to original
SAVE_RELABELLED_COO = "SPANNER/relabelled_coo.npy"  # Path to save the relabelled COO

# Read the input file
edgelist = pd.read_csv(INPUT_PATH, sep="\t", header=None)
edgelist = np.array(edgelist)

# Initialize arrays and variables
inverse = np.zeros(int(edgelist[0, 0]))
mappings = np.ones(int(edgelist[1:].max()) + 1) * -1
sanitized = []
sanitized.append([edgelist[0, 0], edgelist[0, 1]])
counter = 0

# Sanitize the edgelist
for i in tqdm(range(1, len(edgelist))):
    e = edgelist[i]
    u = int(e[0])
    v = int(e[1])
    if mappings[u] == -1:
        mappings[u] = counter
        inverse[counter] = u
        counter += 1
    if mappings[v] == -1:
        mappings[v] = counter
        inverse[counter] = v
        counter += 1
    x = int(mappings[u])
    y = int(mappings[v])
    sanitized.append([x, y])

# Process the sanitized data
coo = np.array(sanitized)
coo_u = np.sort(coo[1:], axis=1)
coo_u = np.unique(coo_u, axis=0)
row_indices = np.where(coo_u[:, 0] == coo_u[:, 1])[0]
if len(row_indices):
    coo_u = np.delete(coo_u, row_indices, axis=0)
edges = len(coo_u)
print(edges, np.max(coo_u) + 1)
vertices = np.max(coo_u) + 1
coo = np.vstack(([vertices, edges], coo_u))

# Save the mappings and relabelled COO to .npy files
np.save(SAVE_MAPPING_TO_OG_COO, inverse)
np.save(SAVE_RELABELLED_COO, coo)

# Save the sanitized edgelist to a text file
np.savetxt(OUTPUT_PATH, coo, fmt='%d', delimiter="\t")

