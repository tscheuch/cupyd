## A collection of hopefully useful thoughts

- The CHD package is used to store the boundary conditions per cell per stress period of the water table
  (altura de agua).

- The DIS package holds the **itmuni** attribute which represents the **Time units, default is days (4)** per stress period.

- In order to change the boundary conditions, we would need to rewrite the CHD package every time a stress period ends.
That can be made following these steps:

  1. Change the **stress_period_data** 0-layer to the correct stress period data.
  2. Call the **write_file** function.
  3. Run the simulation.
