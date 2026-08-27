import torch


def load_gro_file(fname):
    with open(fname, 'r') as f:
        name = f.readline()
        natom = int(f.readline())
        p_mol = torch.zeros(natom)
        p_type = torch.zeros(natom)
        p_position = torch.zeros(natom, 3)
        residual_index = {}
        next_residual_index = 0
        for i in range(natom):
            residual, type, index, x, y, z = f.readline().split()
            if residual not in residual_index.keys():
                residual_index[residual] = next_residual_index
                next_residual_index += 1
            assert type in ["OW", "HW1", "HW2"]

            p_mol[i] = residual_index[residual]
            p_type[i] = 0 if type == "OW" else 1
            p_position[i] = torch.tensor([float(c) for c in (x, y, z)])
        lx, ly, lz = (float(c) for c in f.readline().split())
        assert lx == ly == lz
        assert (natom % 3) == 0
        assert (p_type == 0).type(torch.int).sum() == natom//3
        assert (p_type == 1).type(torch.int).sum() == 2*natom//3
        # TODO: some more checks that right number of atoms in residues
        return p_mol, p_type, p_position, lx


p_mol, p_type, p_position, L = load_gro_file("test.gro")
print(p_mol, p_type, p_position, L)
