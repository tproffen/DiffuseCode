# API reference

```{eval-rst}
.. currentmodule:: diffuse
```

## Structures

```{eval-rst}
.. autofunction:: read_structure

.. autofunction:: make_crystal

.. autofunction:: from_ase

.. autoclass:: Crystal
   :members: atoms, atom_types, positions, cartesian, lattice_vectors, to_ase, view, volume, a, label

.. autoclass:: StructureView
   :members: show, write_html
```

## Powder diffraction

```{eval-rst}
.. autofunction:: powder

.. autoclass:: PowderPattern
   :members: two_theta, plot
```

## Single-crystal diffraction

```{eval-rst}
.. autofunction:: fourier

.. autoclass:: FourierMap
   :members: hkl, plot, is_line
```

## Pair distribution function

```{eval-rst}
.. autofunction:: pdf

.. autoclass:: PDFPattern
   :members: plot
```

## Running DISCUS directly

```{eval-rst}
.. autofunction:: run_macro

.. autofunction:: session

.. autoclass:: Session
   :members: run, value, dataset, dataset_2d, workdir, last_output, last_error,
             reset, seed, random_state, set_random_state, backend
```

## State and random numbers

```{eval-rst}
.. autofunction:: seed

.. autofunction:: reset
```

## Errors

```{eval-rst}
.. autoexception:: DiffuseError

.. autoexception:: MacroError
```
