# Output Files

### Facets Data File

For a typical `LDPMgeo000-data-facets.dat` file:

```
// Data Structure:
// Tet IDx IDy IDz Vol pArea cx cy cz px py pz qx qy qz sx sy sz mF
// One line per facet, ordering is Tet 1 (Facet 1-12),...,Tet N (Facet 1-12)
// Note: All indices are zero-indexed
```

* **Tet:** Index of tetrahedra (_integer_)
* **IDx IDy IDz:** Indices of facet vertices (_ints_)
* **Vol**: Tet volume associated to facet (_float_)
* **pArea**: Projected area of facet (_float_)
* **cx, cy, cz**: Coordinates of facet center (floats)
* **px, py, pz**: Normal to projected facet (_floats_)
* **qx, qy, qz**: First projected tangential direction (_floats_)
* **sx, sy, sz**: Second projected tangential direction (_floats_)
* **mF**: Material Flag (_integer_); where,
  * 0 = Single-Phase LDPM
  * 1 = ITZ (for P-LDPM)
  * 2 = Binder (for P-LDPM)
  * 3 = Aggregate (for P-LDPM)

### Face Facets Data File

For a typical `LDPMgeo000-data-faceFacets.dat` file:

```
// Data Structure:
// n x1 y1 z1 x2 y2 z2 x3 y3 z3
// One line per face facet, first number is the tet node index
// Note: All indices are zero-indexed
```

* **n:** Tet node index (_int_)
* **x1:** X Coordinate for vertex 1 (_float_)
* **y1:** Y Coordinate for vertex 1 (_float_)
* **z1:** Z Coordinate for vertex 1 (_float_)
* **x2:** X Coordinate for vertex 2 (_float_)
* **y2:** Y Coordinate for vertex 2 (_float_)
* **z2:** Z Coordinate for vertex 2 (_float_)
* **x3:** X Coordinate for vertex 3 (_float_)
* **y3:** Y Coordinate for vertex 3 (_float_)
* **z3:** Z Coordinate for vertex 3 (_float_)
