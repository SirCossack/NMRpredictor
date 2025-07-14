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




import mysql.connector as sql

connection = sql.connect(user='root', host='127.0.0.1', raise_on_warnings=True, database='nmrshiftdb')
cursor = connection.cursor()


cursor.execute('WITH Spec AS ( SELECT s.molecule_id, s.spectrum_id, s.spectrum_type_id, nc.value AS condition_value FROM Spectrum_condition sc JOIN nmr_condition nc ON sc.condition_id = nc.condition_id JOIN spectrum s ON s.spectrum_id = sc.spectrum_id WHERE nc.condition_type_id=3),'
               ' Specc AS (SELECT st.name, st.Spectrum_type_id FROM spectrum_type_condition stc JOIN spectrum_type st on stc.spectrum_type_id = st.spectrum_type_id JOIN condition_type ct on stc.condition_type_id=ct.condition_type_id  WHERE ct.condition_type_id=3),'
               ' Signals AS (SELECT sa.atom_id, sh.signal_id, sh.axis, sh.value AS shift_value, ns.spectrum_id, ns.intensity, ns.multiplicity FROM shift sh JOIN nmr_signal ns ON sh.signal_id = ns.signal_id JOIN signal_atom sa ON sh.signal_id = sa.signal_id),'
               ' Bonds AS (SELECT b.bond_id, a.atom_id FROM bond_atom ba JOIN atom a on ba.atom_id = a.atom_id JOIN bond b on ba.bond_id = b.bond_id)'
               ' SELECT Bonds.atom_id, Bonds.bond_id, Molecule.smiles_string, spec.spectrum_id, Signals.signal_id, axis, shift_value, intensity, multiplicity from Spec JOIN Specc ON Spec.Spectrum_type_id = Specc.Spectrum_type_id JOIN Signals ON Signals.spectrum_id = Spec.spectrum_id JOIN Molecule ON Molecule.molecule_id = Spec.molecule_id JOIN Bonds ON Signals.atom_id = Bonds.atom_id WHERE name="1H" ORDER BY Spec.spectrum_id')
result = cursor.fetchall()
cols = [x[0] for x in cursor.description]

id = None
for r in result:
    r = dict(zip(cols, r))

    if id == None:
        id = r['spectrum_id']

    cursor.execute(f'SELECT * from atom where atom_id = {r["atom_id"]}')
    print(cursor.fetchall())
    if r['spectrum_id'] != id:
        break
    print(r)


def input_returner(HOSE: str) -> list:
    """
    :param HOSE: HOSE string
    :return: List of 80 values to be put into neural network
    """

    return []