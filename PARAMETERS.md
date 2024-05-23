### `--landuse`

Land use is primarily used on IFM to estimate how quick water left on the surface will flow.
Such an estimate is given by Manning's roughness coefficients, also called _Manning's N_.
There are [several resources](http://www.lmnoeng.com/manningn.htm) on the literature that
estimate values for common land cover types.

Here we aggregate several estimates for Manning's N and incorporate landuse classes that are
defined by the OpenStreetMap project. When defining a fixed value for `landuse` you may pick
the one that more closely represents the geographical region being modeled.

| OSM _landuse_ class   | Manning's N | Observation           |
| ------------------------|-------------|-----------------------|
| allotments              | 0.15        | Clover                |
| basin                   | 0.0137      | Developed/industrial  |
| brownfield              | 0.02        | Graveled surface      |
| cemetery                | 0.012       | Graveled surface      |
| commercial              | 0.0137      | Developed/industrial  |
| construction            | 0.0137      | Developed/industrial  |
| farm                    | 0.246       | Cotton/Soy            |
| farmland                | 0.246       | Cotton/Soy            |
| farmyard                | 0.246       | Cotton/Soy            |
| forest                  | 0.192       | Forest                |
| garages                 | 0.0137      | Developed/industrial  |
| grass                   | 0.45        | Grass (bluegrass sod) |
| greenfield              | 0.235       | Pasture               |
| greenhouse_horticulture | 0.15        | Clover                |
| industrial              | 0.0137      | Developed/industrial  |
| landfill                | 0.02        | Graveled surface      |
| meadow                  | 0.24        | Dense grass           |
| military                | 0.150       | Sparsely vegetated    |
| orchard                 | 0.235       | Pasture               |
| plant_nursery           | 0.41        | Bermuda grass         |
| quarry                  | 0.320       | Gullied land          |
| railway                 | 0.0137      | Developed/industrial  |
| recreation_ground       | 0.05        | Bare field            |
| reservoir               | 0.035       | [Major rivers](http://www.lmnoeng.com/manningn.htm)|
| residential             | 0.0137      | Developed/industrial  |
| retail                  | 0.0137      | Developed/industrial  |
| salt_pond               | 0.035       | [Major rivers](http://www.lmnoeng.com/manningn.htm)|
| village_green           | 0.40        | Lawns                 |
| vineyard                | 0.09        | Rows crops            |


| OSM _natural_ class | Mannning's N | Observation              |
|---------------------|--------------|--------------------------|
| arete               |    0.035     | Stony Cobbles            |
| bay                 |    0.02      | Bare sand                |
| beach               |    0.02      | Bare sand                |
| cliff               |    0.035     | Stony Cobbles            | 
| coastline           |    0.02      | Bare sand                |
| fell                |    0.17      | Dense grass              |
| grassland           |    0.24      | Dense grass              |
| glacier             |    0.03      | (guessed)                |
| heath               |    0.05      | Bare field - no residue  |
| mud                 |    0.025     | Bare clay-loam           |
| ridge               |    0.15      | Grass and pasture        |
| sand                |    0.02      | Bare sand                |
| scree               |    0.02      | Graveled surface         |
| scrub               |    0.320     | Gullied land             |
| sinkhole            |    0.320     | Gullied land             |
| spring              |    0.035     | Major rivers             |
| stone               |    0.035     | Stony Cobbles            |
| tree                |    0.192     | Forest                   |
| tree_row            |    0.40      | Dense growth             |
| water               |    0.035     | Major rivers             |
| wetland             |    0.035     | Major rivers             |
| wood                |    0.184     | Forest                   |


### `--soil-ep`, `--soil-hc`, `--soil-ph`

The following values have been taken from Rawls, W., Brakensiek, D., and Miller, N. (1983),
"Green-ampt Infilitration Parameters from Soils Data.", J. Hydraul. Eng., 109(1), 62-70.

Some values are not defined in the table below but can be determined as a combination of
similar texture classes.

| Soil texture class | Effective porosity | Hydraulic Conductivity | Pressure head |
|--------------------|--------------------|------------------------|---------------|
| clay (heavy)       | 0.385              | 0.06                   | 31.63         |
| silty clay         | 0.423              | 0.10                   | 29.22         |
| clay (light)       | NA                 | NA                     | NA            |
| silty clay loam    | 0.432              | 0.20                   | 27.30         |
| clay loam          | 0.390              | 0.20                   | 20.88         |
| silt               | NA                 | NA                     | NA            |
| silt loam          | 0.486              | 0.68                   | 16.68         |
| sandy clay         | 0.330              | 0.30                   | 21.85         |
| loam               | 0.434              | 1.32                   | 8.89          |
| sandy clay loam    | 0.330              | 0.30                   | 21.85         |
| sandy loam         | 0.412              | 1.15                   | 11.01         |
| loamy sand         | 0.401              | 1.30                   | 6.13          |
| sand               | 0.417              | 1.45                   | 4.95          |

Some important notes:
1. Effective porosity is given in cm3.
2. Hydraulic conductivity at saturation is given in cm/h. For IFM's Green-Ampt soil model,
   please divide the value seen in the table by 2.
3. Wetting front suction head (capillary head) is given in cm.

