import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import dijkstra
from ase.neighborlist import NeighborList
from ase.data import covalent_radii
from ase.symbols import symbols2numbers



def check_ring_is_periodic(ring, offsets):
    ''' 
    Check if the ring wraps around the period cell, i.e., is not a true ring.
    Args:
        ring (list(int)): Atom indices of the ring
        offsets (list(np.ndarray)): Unit cell offsets for all atoms
    Returns:
        bool: True if the ring is periodic, False otherwise
    '''
    total_offset = np.zeros(3)
    for i in range(len(ring) - 1):
        total_offset += offsets[ring[i], ring[i+1]]
    total_offset += offsets[ring[-1], ring[0]]
    return np.all(total_offset == 0)

def find_rings(ats, radii_factor=1.3, repeat=(1, 1, 1), bonds=None):
    ''' Find rings in the unit cell.
        Rings are found according to the definition of L. Guttman, J. Non-Cryst. Solids 1990, 116.
        Args:
            ats (ase.Atoms): Atoms object containing the structure
            radii_factor (float): Factor to multiply covalent radii for neighbor search
            repeat (tuple(int, int, int)): How often to repeat the unit cell in each direction. Increase for small cells.
            bonds (list(tuple(str, str))): List of allowed bonds, e.g., [('C', 'C'), ('C', 'O')], can be None to allow all bonds.
    '''
    s = ats.repeat(repeat)
    pos = s.get_positions()
    nat = len(s)
    lat = s.get_cell()
    els = s.get_chemical_symbols()
    radii = covalent_radii[symbols2numbers(els)]

    if bonds is not None:
        # Don't need to find neighbors for elements not included in bonds
        elements = set().union(*bonds)
        radii = [x if el in elements else 0. for el, x in zip(els, radii)]
        radii = np.array(radii, dtype=float)
    
    nl = NeighborList(radii * radii_factor, self_interaction=False, bothways=False, skin=0.)
    nl.update(s)

    # pairwise distances
    d = np.zeros((nat, nat))
    # unit cell offsets for all atoms
    all_offsets = np.zeros((nat, nat, 3), dtype=int)

    for i in range(nat):
        indices, offsets = nl.get_neighbors(i)
        rs = pos[indices, :] + offsets @ lat - pos[i, :]
        ds = np.linalg.norm(rs, axis=1)
        d[i, indices] = ds
        d[indices, i] = ds
        all_offsets[i, indices] = offsets
        all_offsets[indices, i] = -offsets

        # set all neighbors that are not allowed to 0 (will be removed later)
        if bonds is not None:
            for j in indices:
                if not ((els[i], els[j]) in bonds or (els[j], els[i]) in bonds):
                    d[i, j] = 0
                    d[j, i] = 0

    # sparse matrix of bonds, removes zero entries
    d = csr_array(d)

    # now find the rings 
    rings = {}
    for i in range(len(ats)):
        indices, offsets = nl.get_neighbors(i)
        for j, _ in zip(indices, offsets):
            if bonds is not None: 
                # skip if bond is not allowed
                if not ((els[i], els[j]) in bonds or (els[j], els[i]) in bonds):
                    continue
            # Remove the bond between the two selected neighbors
            d_tmp = d.copy()
            d_tmp[i, j] = 0
            d_tmp[j, i] = 0
            d_tmp.eliminate_zeros()
            # Now find the shortest path between i and j
            dist_matrix, predecessors = dijkstra(d_tmp, indices=i, return_predecessors=True, directed=False, unweighted=True, limit=np.inf)
            if dist_matrix[j] < np.inf:
                k = j
                ring = [k]
                while predecessors[k] != i:
                    k = predecessors[k]
                    ring.append(k)
                ring.append(i)

                if not check_ring_is_periodic(ring, all_offsets):
                    print('WARNING: ring is wrapping around periodic cell! Consider increasing `repeat`.')
                    continue

                ring = [x % len(ats) for x in ring]  # take it back to primary cell
                rings[tuple(sorted(ring))] = ring


    return list(rings.values())
