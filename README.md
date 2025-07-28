# Rings

This is a little Python script to find rings in atomistic structures.

It is using the ring definition according to [Guttman, L. (1990). Ring structure of the crystalline and amorphous forms of silicon dioxide. Journal of Non-Crystalline Solids, 116(2), 145-147](https://doi.org/10.1016/0022-3093(90)90686-G).

The function should be straight forward to use. An example is given in `main.py`.

```python
    from rings import find_rings
    from ase.io import read, write

    # Read input structure 
    ats = read('test.extxyz')

    # Find the rings
    rings = find_rings(ats, radii_factor=1.3, bonds=[('P', 'O'), ('Al', 'O')], repeat=(2, 2, 2))
    print(f'Found {len(rings)} rings')

    # `rings` now contains a list of index lists, identifying the atoms participating in each ring

    # Create a list of structures containing only the rings
    rs = []
    for r in rings:
        ss = ats[r]
        rs.append(ss)

    # Save it to have a look at it
    write('rings.extxyz', rs)
```



