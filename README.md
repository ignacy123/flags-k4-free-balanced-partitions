This repository contains files for computer-assisted parts of the proof for the paper TODO.

# Requirements
To run calculations in this repository, you need:
- `csdp` -- compiled with more lenient exit conditions. You can find a repository here: TODO.
- `sage`
- `cmake`.

# Setup
```
./init-cmake.sh
```

This will setup files necessary to run calculations and prepare metadata that will allow LSP completions to work.

# Running SDPs
```
cmake --build build --target <problem-name> --config Release
```
Problem names can be found in `CMakeLists.txt`.

# Rounding
TODO. Files are in the repository, commands have to be added...
