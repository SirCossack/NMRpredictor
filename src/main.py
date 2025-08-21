"""
First commit probably, as i escaped SQL hell. Shouldn't need anything more.
The plan is to have a Pytorch Neural Network and input 80 values into it, expecting 1 NMR shift value in return.
80 values is 40 atoms in total, assuming:
 - only 1 atom directly connected to H (2 inputs, one for atom type and bond)
 - max 3 atoms connected to H via 2 bonds (6 inputs, 3 atom types and 3 bonds connecting)
 - max 9 atoms connected via 3 bonds (18 inputs)
 - max 27 atoms connected via 4 bonds (54 inputs)
Atoms connected further will be omitted.

Now I gotta figure out PyTorch's neural network input scheme and make a function to transform HOSE code into these 80 inputs
"""


import torch.nn as nn
import mysql.connector as sql
import torch
import torch.optim as op
from torch.utils.data import Dataset, DataLoader
from rdkit import Chem
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdmolops import AddHs


class Network(nn.Module):
    def __init__(self):
        super().__init__()
        self.lay1 = nn.Linear(13*121, 100)
        self.activation = nn.ReLU()
        self.lay2 = nn.Linear(100, 1)

    def forward(self, x):
        x = self.lay1(x)
        x = self.activation(x)
        x = self.lay2(x)
        x = self.activation(x)
        return x

class NMRData(Dataset):
    def __init__(self, data, targets):
        if data.size()[0] != len(targets):
            raise Exception("Uneven number of data tensors and targets")
        self.data = data
        self.targets = targets

    def __len__(self):
        return len(self.targets)

    def __getitem__(self, idx):
        return (self.data[idx].float(), self.targets[idx])

def HOSE_to_list(HOSE:str, smiles) -> list:
    layer = 0
    layerdict = {0: 1,
                 1: 4,
                 2: 13,
                 3: 40,
                 4: 121}

    if "H" not in HOSE:
        pass


    bond = "-"
    tensor = []
    bonddict = {"-": [1,0,0,0],
                '*': [0,1,0,0],
                '=': [0,0,1,0],
                "%": [0,0,0,1]}
    atomdict = {"C": [1,0,0,0,0,0,0,0,0],
                "H": [0,1,0,0,0,0,0,0,0],
                "O": [0,0,0,0,0,0,0,1,0],
                "N": [0,0,1,0,0,0,0,0,0],
                "X": [0,0,0,1,0,0,0,0,0],
                "Y": [0,0,0,0,1,0,0,0,0],
                "I": [0,0,0,0,0,1,0,0,0],
                "S": [0,0,0,0,0,0,1,0,0],
                "F": [0,0,0,0,0,0,0,0,1]}

    for char in HOSE:
        if char in {'-', '=', '*', "%"}:
            bond = char

        elif char == ",":
            pass

        elif char.isalpha() or char == "&":
            atomtensor = atomdict.get(char)
            if atomtensor:
                tensor.extend(atomtensor)
            else:
                tensor.extend(list([0 for x in range(len(atomdict["C"]))]))
            tensor.extend(bonddict[bond])
            bond = "-"

        elif char in {"/", "(", ")"}:
            tensor.extend(list([0 for x in range(layerdict[layer] * (len(bonddict["-"]) + len(atomdict["C"])) - len(tensor))]))
            if layer == 4:
                #tensor = torch.tensor(tensor)
                #tensor = tensor.view(layerdict[layer], 13)
                return tensor
            layer += 1

        elif char in {"+", "-", "@", "#", "|", "\\"}:
            pass

        else:
            print("unknown character:, ",char)
            print("hose code: ", HOSE)
            print("moleucle: ", smiles)

    layer = 4
    tensor.extend(list([0 for x in range(layerdict[layer] * (len(bonddict["-"]) + len(atomdict["C"])) - len(tensor))]))
    #tensor = torch.tensor(tensor)
    #tensor = tensor.view(layerdict[layer], 13)
    return tensor

def list_to_tensor(tensor) -> torch.Tensor:
    tensor = torch.tensor(tensor)
    tensor = tensor.view(121, 13)
    return tensor

def Atom_to_HOSE(_atoms, _seen=None, layer=0 ):
    if not _seen:
        _seen = []
    if layer == 6:
        return ""
    if type(_atoms) != list:
        _atoms = [_atoms]

    HOSE = ""

    bonddict = {BondType.AROMATIC: "*",
                BondType.SINGLE: "",
                BondType.DOUBLE: "=",
                BondType.TRIPLE: "%"}

    atomdict = {"Cl": "X",
                 "Br": "Y",
                 "Si": "Q",
                 "C": "C",
                 "P": "P",
                 "H": "H",
                 "O": "O",
                 "F": "F",
                 "N": "N",
                 "S": "S",
                 "I": "I"}
    atoms = []
    for atom in _atoms:
        for bond in atom.GetBonds():
            otheratom = bond.GetOtherAtom(atom)
            if otheratom.GetIdx() in _seen:
                continue
            _seen.append(otheratom.GetIdx())
            bondtype = bond.GetBondType()
            if atomdict.get(otheratom.GetSymbol()):
                HOSE += bonddict[bondtype] + atomdict[otheratom.GetSymbol()]
            else:
                HOSE += bonddict[bondtype] + "&"
            atoms.append(otheratom)
        HOSE += ","
    HOSE = HOSE[:-1]
    HOSE += "/"

    return HOSE + Atom_to_HOSE(atoms, _seen, layer+1)


### Query the DB
connection = sql.connect(user='root', host='127.0.0.1', raise_on_warnings=True, database='nmrshiftdb')
cursor = connection.cursor()
cursor.execute('WITH Spec AS ( SELECT s.molecule_id, s.spectrum_id, s.spectrum_type_id, nc.value AS condition_value FROM Spectrum_condition sc JOIN nmr_condition nc ON sc.condition_id = nc.condition_id JOIN spectrum s ON s.spectrum_id = sc.spectrum_id WHERE nc.condition_type_id=3),'
               ' Specc AS (SELECT st.name, st.Spectrum_type_id FROM spectrum_type_condition stc JOIN spectrum_type st on stc.spectrum_type_id = st.spectrum_type_id JOIN condition_type ct on stc.condition_type_id=ct.condition_type_id  WHERE ct.condition_type_id=3),'
               ' Signals AS (SELECT sa.atom_id, sh.signal_id, sh.axis, sh.value AS shift_value, ns.spectrum_id, ns.intensity, ns.multiplicity FROM shift sh JOIN nmr_signal ns ON sh.signal_id = ns.signal_id JOIN signal_atom sa ON sh.signal_id = sa.signal_id),'
               ' Bonds AS (SELECT a.hose_code_with_rings, b.bond_id, a.atom_id FROM bond_atom ba JOIN atom a on ba.atom_id = a.atom_id JOIN bond b on ba.bond_id = b.bond_id)'
               ' SELECT Bonds.hose_code_with_rings, Bonds.atom_id, Bonds.bond_id, Molecule.smiles_string, spec.spectrum_id, Signals.signal_id, axis, shift_value, intensity, multiplicity from Spec JOIN Specc ON Spec.Spectrum_type_id = Specc.Spectrum_type_id JOIN Signals ON Signals.spectrum_id = Spec.spectrum_id JOIN Molecule ON Molecule.molecule_id = Spec.molecule_id JOIN Bonds ON Signals.atom_id = Bonds.atom_id WHERE name="1H" ORDER BY Spec.spectrum_id')
result = cursor.fetchall()

### i should probably cache the sql query

### Try to use gpu for fancy shmancy neuromancy --- maybe later
print("Finished database query...")

### Construct the Dataset
cols = [x[0] for x in cursor.description]
data = []
targets = []
molecule_tensors = {}
for r in result:
    r = dict(zip(cols, r))
    smiles = str(r['smiles_string'])[2:-1]
    """this is for input validation but doesnt work since HOSE in database is in slightly different order than HOSE from rdkit.chem molecules and my function
    molecule = Chem.MolFromSmiles(smiles)
    if not molecule:
        continue
    molecule = AddHs(molecule)
    if molecule not in molecule_tensors:
        molecule_tensors[molecule] = []
        for atom in molecule.GetAtoms():
            if atom.GetSymbol() == "H":
                hose = Atom_to_HOSE(atom)
                _tensor = HOSE_to_list(hose, smiles)
                molecule_tensors[molecule].append(_tensor)
    """

    hose = r['hose_code_with_rings'][r['hose_code_with_rings'].find(";")+1:]
    tensor = HOSE_to_list(hose, smiles)

    """this is for input validation but doesnt work since HOSE in database is in slightly different order than HOSE from rdkit.chem molecules and my function
    if tensor in molecule_tensors[molecule]:
        tensor = list_to_tensor(tensor)
        data.append(tensor)
        shift = r['shift_value']
        targets.append(shift)
    else:
        pass #print(tensor, molecule_tensors[molecule], tensor in molecule_tensors[molecule])
        """
    tensor = list_to_tensor(tensor)
    data.append(tensor)
    shift = r['shift_value']
    targets.append(shift)

print('finished constructing dataset')

targets = torch.tensor(targets).float()
data = torch.stack(data)
nmr_data = NMRData(data, targets)
batchsz = 25
training_data, eval_data, test_data = torch.utils.data.random_split(nmr_data, (0.8, 0.05, 0.15))
training_data = DataLoader(training_data, batch_size=batchsz, shuffle=True)
test_data = DataLoader(test_data, batch_size=batchsz, shuffle=True)
eval_data = DataLoader(eval_data, batch_size=batchsz, shuffle=True)

### Initialize and train the network
network = Network()
epochs = 10
loss_func = nn.MSELoss()
learning_rate = 0.01
optimizer = op.SGD(network.parameters(), lr=learning_rate)



for epoch in range(epochs):
    network.train()
    for data, target in training_data:
        data = data.reshape(data.shape[0], -1)
        outputs = network(data)
        outputs = outputs.squeeze()
        loss = loss_func(outputs, target)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    network.eval()
    ppm1 = 0
    ppm05 = 0
    ppm01 = 0
    all = 0
    for data, target in eval_data:
        data = data.reshape(data.shape[0], -1)
        outputs = network(data)
        outputs = outputs.squeeze()
        for i in range(len(outputs)):
            all +=1
            diff = abs(outputs[i].item() - target[i].item())
            if diff < 0.1:
                ppm01 +=1
            elif diff < 0.5:
                ppm05 += 1
            elif diff < 1:
                ppm1 +=1
    print(f"Epoch {epoch}: {ppm01/all:.2f} correct-ish, {ppm05/all:.2f} within 0,5ppm and {ppm1/all:.2f} within 1ppm. ")
