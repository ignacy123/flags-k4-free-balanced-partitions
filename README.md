This repository contains files for computer-assisted parts of the proof for the paper TODO.

A large part of this codebase is originally based on Bernard Lidický's work, but it was modified and developed by the authors of this paper.

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

Problem names can be found in `CMakeLists.txt`.
# Running SDPs
```
cmake --build build --target <problem-name>-solve-approximate --config Release
```

# Rounding
```
cmake --build build --target <problem-name>-round --config Release
```
