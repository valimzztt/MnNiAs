import numpy as np
import json
from pymatgen.core.structure import Structure
from smol.io import load_work
from smol.cofe import ClusterSubspace, StructureWrangler, ClusterExpansion, RegressionData
from smol.io import save_work


file_path = 'MnNiAs-smol/ce-data/ce_MnNiAs_smol.mson'

work = load_work(file_path)
for name, obj in work.items():
    print(f'{name}: {type(obj)}\n')
expansion  = work.get("ClusterExpansion")
wrangler = work.get("StructureWrangler")
print(expansion)
