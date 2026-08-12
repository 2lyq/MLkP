# Regional_MSRs

Regional_MSRs is a C++ implementation for mining regional co-location / MSR patterns from spatial feature instances.

## Requirements

- CMake 4.0 or newer
- A C++20 compiler

The project has been tested with the MinGW toolchain bundled with CLion on Windows. On Linux or macOS, the generated executable name will usually be `Regional_MSRs` instead of `Regional_MSRs.exe`.

## Input CSV format

The input file is read as a CSV file with a header row. Each record is expected to contain at least four columns:

```csv
Feature,InstanceID,LocationX,LocationY
A,1,116.123456,39.123456
B,1,116.223456,39.223456
```

Column meaning:

- `Feature`: feature type or category
- `InstanceID`: instance id within the feature
- `LocationX`: x coordinate
- `LocationY`: y coordinate

## Build

From the project root:

```bash
cmake -S . -B cmake-build-release -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build-release
```

If you use CLion, open this folder as a CMake project and build the `Regional_MSRs` target.

## Run

```bash
./cmake-build-release/Regional_MSRs.exe <dataset path> <minimum prevalent index> <candidate ratio> <distance> <k>
```

Example:

```bash
./cmake-build-release/Regional_MSRs.exe data/sample.csv 0.4 0.02 500 3
```

Parameters:

- `dataset path`: input CSV file path
- `minimum prevalent index`: minimum prevalence threshold for global co-location filtering
- `candidate ratio`: candidate filtering ratio
- `distance`: spatial neighbor distance threshold
- `k`: pattern size

## Output files

The program writes result files to the current working directory:

- `PatternSummary<k>.txt`: k-size local co-location pattern count, k-size global co-location pattern count, and runtime
- `Colocation.txt`: global co-location patterns
- `PatternsScore<k>.txt`: regional pattern scores
- `Regional Colocation_/<k>/*.csv`: mined regional instances for each pattern/MSR

