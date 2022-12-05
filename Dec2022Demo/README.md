# Demo for December 2022

The organization of this folder is the following:

- `outputs` for json created by scripts
- `scripts` for scripts meant to test out functionality
- `wiring_diagrams` for scripts that define wiring diagrams
- `notebooks` for Pluto or Jupyter notebooks meant for presenting

What should *not* be in this folder?

- Almost all function definitions should be in `../lib`. If it's worth abstracting into a function, it's worth doing it right and putting it in a library. Function definitions in `../lib` should be grouped sensibly into modules.
- Basic models like `SIR`, `SIRD`, etc., and their typings. These should also be in `../lib`.
