"""Isolated pymatgen adsorption-site enumeration worker.

This module is executed in a subprocess so that pymatgen/ruamel import crashes
cannot take down the main CatDT workflow process.
"""

from __future__ import annotations

import json
import sys

from ase.io import read
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core import Lattice
from pymatgen.io.ase import AseAtomsAdaptor


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: pymatgen_site_worker.py '<payload-json>'")

    payload = json.loads(sys.argv[1])
    slab = read(payload["slab_path"])
    anchor = payload["anchor_pos"]
    min_dist = float(payload["min_dist"])
    max_dist = float(payload["max_dist"])

    pmg_slab = AseAtomsAdaptor.get_structure(slab)
    asf = AdsorbateSiteFinder(pmg_slab)
    raw_sites = asf.find_adsorption_sites(symm_reduce=0)
    lattice = Lattice(pmg_slab.lattice.matrix)

    def pbc_distance(cart_a, cart_b):
        frac_a = lattice.get_fractional_coords(cart_a)
        frac_b = lattice.get_fractional_coords(cart_b)
        return float(lattice.get_distance_and_image(frac_a, frac_b)[0])

    sites = []
    seen_xy = set()
    for stype in ("hollow", "bridge", "ontop"):
        for pos3d in raw_sites.get(stype, []):
            dist = pbc_distance(pos3d, anchor)
            if not (min_dist <= dist <= max_dist):
                continue
            key = (round(float(pos3d[0]), 1), round(float(pos3d[1]), 1))
            if key in seen_xy:
                continue
            seen_xy.add(key)
            sites.append(
                {
                    "position": [float(pos3d[0]), float(pos3d[1]), float(pos3d[2])],
                    "dist_from_anchor": dist,
                    "site_label": stype,
                }
            )

    if not sites:
        for stype in ("hollow", "bridge", "ontop"):
            for pos3d in raw_sites.get(stype, []):
                dist = pbc_distance(pos3d, anchor)
                if 2.0 <= dist <= 6.0:
                    sites.append(
                        {
                            "position": [float(pos3d[0]), float(pos3d[1]), float(pos3d[2])],
                            "dist_from_anchor": dist,
                            "site_label": stype,
                        }
                    )

    sites.sort(key=lambda s: float(s.get("dist_from_anchor", 1.0e9)))
    print(json.dumps(sites))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
