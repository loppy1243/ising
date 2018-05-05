# Very basic 2D Ising model simulation

This is a git repository. The current branch, if you've just decompressed this directory,
should be called `factoring`. This is the nice, pretty, commented, and documented version of
the code.

However, it runs like molasses (comparatively).

The `master` contains the fast version of the code, which I don't want to touch. If necessary,
you can probably extrapolate the documentation from the nice code to the bad code.

The directory `sims` contains all the simulations I generated, along with a file "INDEX.txt"
that was for note-jotting as I went along. This is likely illegible. The nicer file is
"PRESENT.txt", which is the list I made for what I planned to show in my presentation. It is,
at least, legible.

To run the program do
```
% julia Main.jl
```
if on "factoring" or
```
% julia IsingMain.jl
```
if on `master`. Modifying the main() function in `Main.jl` (`IsingMain.jl`) to run different
simulations should be fairly straightforward.
