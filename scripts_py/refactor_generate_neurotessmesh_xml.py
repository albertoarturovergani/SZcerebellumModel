from bsb.core import from_storage
import numpy as np
import xml.etree.ElementTree as ET

# read network and cell types
network = from_storage("mouse_cereb_nest.hdf5")

# Set output file name

outname = "NTM_out.xml"

ct = network.cell_types

cell_type = []
cell_layer = []
cell_pos = []
cell_morpho = []
cell_morpho_ids = []

# cell id
gid = 0

# layer id (= population id for neurotessmesh)
lid = 0


# Define a ratio for pruning cell numbers in the representation
pruning_factor = dict(
    granule_cell=1000,
    basket_cell=10,
    purkinje_cell=5,
    stellate_cell=10,
    golgi_cell=10,
)
excluded_cell_types = ["glomerulus", "mossy_fibers"]
placement_set_list = [
    pl_set for pl_set in network.get_placement_sets() if pl_set.tag not in excluded_cell_types
]

for ps in placement_set_list:

    cell_morpho_ids.append([])
    print("numb of ", ps.tag, ": ", len(ps))
    # read cells positions
    pos = ps.load_positions()
    # read morphology path
    morpho_swc = (
        ct[ps.tag].get_morphologies()[0].get_meta()["source"].replace(".swc", "_notags.swc")
    )
    cell_morpho.append(morpho_swc)
    added = 0

    for i, p in enumerate(pos):
        subsampleid = int(i / pruning_factor[ps.tag])
        if i % pruning_factor[ps.tag] == 0:
            cell_pos.append(p)
            cell_layer.append(lid)
            cell_morpho_ids[lid].append(added + gid)
            cell_type.append(ps.tag)
            added += 1

    # upd global cell id
    gid += added
    # upd layer id
    lid += 1


# write xml
root = ET.Element("scene", version="0.1")
morphology_se = ET.SubElement(root, "morphology")
columns_se = ET.SubElement(morphology_se, "columns")
neuronmorphologies_se = ET.SubElement(morphology_se, "neuronmorphologies")

# write cell placement
column_se = ET.SubElement(columns_se, "column", id=str("0"))
minicolumn_se = ET.SubElement(column_se, "minicolumn", id=str("0"))

for idx, (layer, pos, ct) in enumerate(zip(cell_layer, cell_pos, cell_type)):
    # geometric tranformation to apply to morphology
    formatted_matrix = (
        "1.0, 0.0, 0.0, "
        + str(pos[0])
        + ", "
        + "0.0, 0.0, 1.0, "
        + str(pos[1])
        + ", 0.0, 1.0, 0.0, "
        + str(pos[2])
        + ", 0.0, 0.0, 0.0, 1.0 "
    )
    transform_se = ET.SubElement(
        minicolumn_se, "neuron", gid=str(idx), layer=str(layer + 1), type=ct
    )
    sub_transform = ET.SubElement(transform_se, "transform").text = formatted_matrix

# write cell id - morphology associations
for i in np.unique(cell_layer):
    # Add morphologies
    str_idx = str(cell_morpho_ids[i]).replace("[", "").replace("]", "").replace(" ", "")
    ET.SubElement(
        neuronmorphologies_se, "neuronmorphology", neurons=str(str_idx), swc=cell_morpho[i]
    )
    print(f"Used {cell_morpho[i]} morphology")

ET.indent(root, "  ")
tree = ET.ElementTree(root)
tree.write(outname)
