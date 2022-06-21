<p align="center">
  <img width="100" height="100" src="https://emojipedia-us.s3.dualstack.us-west-1.amazonaws.com/thumbs/320/apple/325/heart-with-arrow_1f498.png" alt='cupyd'>
</p>

<p align="center"><strong>cupyd --</strong> a romantic integration scheme between SWMM & MODFLOW</p>

---

The integration scheme is composed of 3 modules:
* Spatial Integration of SWMM and MODFLOW elements
* Coupled SWMM and MODFLOW Simulation: Loop for temporal and spatial data exchange
* Post processing and Results Analysis
## Spatial Integration of SWMM and MODFLOW elements
### Inputs:
- SWMM model (swmm_model.inp) **or** Simulation from PySWMM
- SWMM subcatchments **or/and** storage units (S_polygon.shp **or/and** SU_polygon.shp) and junctions (J_points.shp) shapefiles
- MODLFOW model **or** MODFLOW grid + DEM (only for spatial integration)
- List of groundwater zones and junctiones association (optional) 
### Output: 
- gdf with MODFLOW cells associations (elev, S, SU, DRN, drn_to)
- plots
## Coupled SWMM and MODFLOW simulation: Loop for temporal and spatial data exchange 
### Inputs:
- Simulation from pySWMM
- MODFLOW flopy model (ml)
- reporting time steps 
### Outputs:
- Results Time Series
- Sim results (PySWMM)
- Zone Budjet (ZB) results (flopy)
## Post processing and Results analysis
- Continuity Analysis
- SWMM Results 
- MODFLOW Results
- Integration Results




# Library Usage :construction: (for us now)

For the setup and many other things, we have defined many utilities on the Makefile to facilitate your work.
None of them are quite complex, so you are welcome to go and check them or even run the commands separately.

To facilitate the instructions, we use them all around.

1. Ensure that you have python 3.8 installed
2. Install poetry using our custom utilities from the Makefile
```bash
make poetry
```
3. Create a virtual environment and install the project dependencies on it
```bash
make venv-with-dependencies
```


## <a name="commit"></a> Commit Message Format

*This specification is inspired by the [AngularJS commit message format](https://github.com/angular/angular/blob/master/CONTRIBUTING.md#-commit-message-format).*

We will only leave here what we think that is most relevant to our dependency [semantic-release](https://github.com/semantic-release/semantic-release).
#### <a name="commit-header"></a>Commit Message Header

```
<type>: <short summary>
  │            │
  │            └─⫸ Summary in present tense. Not capitalized. No period at the end.              
  │
  └─⫸ Commit Type: build|ci|docs|feat|fix|perf|refactor|test
```

##### Type

Must be one of the following:

* **build**: Changes that affect the build system or external dependencies (example scopes: gulp, broccoli, npm)
* **ci**: Changes to our CI configuration files and scripts (examples: CircleCi, SauceLabs)
* **docs**: Documentation only changes
* **feat**: A new feature
* **fix**: A bug fix
* **perf**: A code change that improves performance
* **refactor**: A code change that neither fixes a bug nor adds a feature
* **test**: Adding missing tests or correcting existing tests

**Please check the commit message implications for the release [here](https://github.com/semantic-release/semantic-release#how-does-it-work)**.
