# API reference

```{eval-rst}
.. currentmodule:: diffuse
```

## Structures

```{eval-rst}
.. autofunction:: read_structure

.. autoclass:: Crystal
   :members: atoms, positions, volume, a
```

## Powder diffraction

```{eval-rst}
.. autofunction:: powder

.. autoclass:: PowderPattern
   :members: two_theta, plot
```

## Running DISCUS directly

```{eval-rst}
.. autofunction:: run_macro

.. autofunction:: session

.. autoclass:: Session
   :members: run, value, dataset, dataset_2d, workdir, last_output
```

## Errors

```{eval-rst}
.. autoexception:: DiffuseError

.. autoexception:: MacroError
```
