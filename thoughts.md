The CHD package is used to store the boundary conditions per cell per stress period of the water table (altura de agua).

The DIS package holds the `itmuni` attribute which represents the `Time units, default is days (4)` per stress period.

We would need to rewrite the CHD package every time a stress period ends in order to change the boundary conditions. That can be made following this steps.

- Change the `stress_period_data` 0 layer to the correct stress period data
- call `write_file`
- Run simulation