# Unscented Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program

Neil Maude, April 2017

---

## Dependencies
Originally provided as SDCND starter code to build under CMake.
Project further developed using Visual Studio 2017 (as a CMake project).
* cmake >= 3.5

## Basic Build Instructions

1. Clone this repo.
2. Compile using your preferred build system 
3. Run executable: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt output.txt`

## Development Notes
Kept the UpdateLidar() and UpdateRadar() functions separate to make the code more 
understandable when referenced in the future (there is scope to consolidate much of the
code within these functions, as the only key difference is the measurement space 
dimension - 2 for lidar, 3 for radar).

## Project Results

### Data file RMSE results
Results when using file `obj_pose-laser-radar-synthetic-input.txt`:

RMSE values (to 3 decimals):
- Px: 0.068 (required < 0.09)
- Py: 0.083 (required < 0.10)
- Vx: 0.280 (required < 0.40)
- Vy: 0.216 (required < 0.30)

### Additional analysis
Additional project analysis can be found in `Unscented Kalman Filter Notes.pdf`.  
This PDF includes:
1. Approximation of `std_a_` and `std_yawdd_` noise parameters.
2. Filter consistency checks using Normalized Innovation Squared (NIS) statistics
3. Performance comparisons between UKF and EKF projects
4. Comparison of lidar-only, radar-only and lidar+radar fusion


