# Very basic 2D Ising model simulation

The `master` branch is the nice, pretty, commented, and documented version of the code and
works with Julia 0.6.

However, it runs like molasses (comparatively).

The `messy` branch contains the fast version of the code (for Julia 0.5 if memory serves),
which I don't want to touch. If necessary, you can probably extrapolate the documentation from
the nice code to the bad code.

A sample output video `sim.mp4` is included.

To run the program do
```
% julia Main.jl
```
if on `master` or
```
% julia IsingMain.jl
```
if on `messy`. Modifying the main() function in `Main.jl` (`IsingMain.jl`) to run different
simulations should be fairly straightforward.
