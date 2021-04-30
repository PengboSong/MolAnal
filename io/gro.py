# coding=utf-8

from gmx.other.input_func import convert_input_type


def readline_gro(line, row):
    assert isinstance(row, int), "Invalid row number."
    assert (isinstance(line, str) and len(line) >=
            44), "Can not parse the given line at row %d." % row
    # Read by column number
    # Molecular id: 1-5, character, right
    mol_id = convert_input_type(line[0:5].strip())
    assert isinstance(
        mol_id, int), "Molecular id at row %d column 1-5 is not an integer." % row

    # Molecular name: 6-8, character, left
    mol_name = convert_input_type(line[5:8].strip())
    assert isinstance(
        mol_name, str), "Molecular name at row %d column 6-8 is not characters." % row

    # Chain identifier: 10, character
    chain = convert_input_type(line[9:10].strip())
    if chain:
        assert isinstance(
            chain, str), "Chain identifier at row %d column 10 is not a character." % row

    # Atom name: 12-15, character, right
    atom_name = convert_input_type(line[11:15].strip())
    assert isinstance(
        atom_name, str), "Atom name at row %d column 12-15 is not characters." % row

    # Atom id: 16-20, integer, right
    atom_id = convert_input_type(line[15:20].strip())
    assert isinstance(
        atom_id, int), "Atom id at row %d column 16-20 is not an integer." % row

    # Coordinate X(A): 21-28, real(8:3), right
    x = convert_input_type(line[20:28].strip())
    assert isinstance(
        x, float), "Coordinate X at row %d column 21-28 is not a real number." % row
    # Coordinate Y(A): 29-36, real(8:3), right
    y = convert_input_type(line[28:36].strip())
    assert isinstance(
        y, float), "Coordinate Y at row %d column 29-36 is not a real number." % row
    # Coordinate Z(A): 37-44, real(8:3), right
    z = convert_input_type(line[36:44].strip())
    assert isinstance(
        z, float), "Coordinate Z at row %d column 37-44 is not a real number." % row

    if len(line.strip()) > 44:
        # Velocity Vx: 45-52, read(8:4), right
        vx = convert_input_type(line[44:52].strip())
        assert isinstance(
            vx, float), "Velocity Vx at row %d column 45-52 is not a real number." % row
        # Velocity Vy: 53-60, read(8:4), right
        vy = convert_input_type(line[52:60].strip())
        assert isinstance(
            vy, float), "Velocity Vy at row %d column 53-60 is not a real number." % row
        # Velocity Vz: 61-68, read(8:4), right
        vz = convert_input_type(line[60:68].strip())
        assert isinstance(
            vz, float), "Velocity Vz at row %d column 61-68 is not a real number." % row
    else:
        vx, vy, vz = 0., 0., 0.

    return {
        "mol_id": mol_id,
        "mol_name": mol_name,
        "chain": chain,
        "atom_name": atom_name,
        "atom_id": atom_id,
        "coordinate": [x, y, z],
        "velocity": [vx, vy, vz]
    }
