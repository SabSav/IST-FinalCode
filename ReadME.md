This project utilizes the Euler method to simulate a 1D Vertex Model, aiming to replicate cell thickness modulation for a wavy cell monolayer, round and cylindric organoids.

## Folder structure:

```
Organoids



    Lutolf: simulations for the geometry u sed in Lutolf paper
        Deterministic:
            DifferentDensity: simulations for different cell densities
            DiffetentShape: simulations for different geometrical parameters
                DifferentLength
                DifferentRadius

        ## Plots Apical Length function of position (z axis):

            ```
            python PlotsALz.py
            ```

            ## Plots Initial and Final Configurations:

            ```
            python PlotsClinderTime.py
            ```

            ## Plots Internuclear Distance (Lutolf paper graph ):

            ```
            python InternuclearDistance.py
            ```
            ## Plots Ratio Thickness Ratio and Modulation of flat and curved region with apical Tension:

            ```
            python ThicknessRatioApicalTension.py
            ```

            ## Plots Modulation change with Density:

            ```
            python ModulationDensity.py
            ```

            ## Plots Modulation change with Length:

            ```
            python ModulationLength.py
            ```

            ## Plots Modulation change with Radius:

            ```
            python ModulationRadius.py
            ```

        Noise:
            DifferentDensity: simulations for different cell densities
            DifferentShape: simulations for different geometrical parameters

            ## Plot Visual Hetorogeneity and YAP:

            ```
            python PlotYAP-Het-Visual.py
            ```

            ## Plot Hetorogeneity in funcion of relative postizion (z axis/L), showing deterministic function and noise contribution:

            ```
            python PlotHet.py
            ```

            ## Plot Apical length in function of relative postizion (z axis/L):

            ```
            python PlotALz.py
            ```

            ## Plot Apical length in function of absolute postizion (z) for different organoid lengths:

            ```
            python PlotALAbsz.py
            ```

            ## Plot Apical length in function of absolute relative postizion (z axis/L) and positional error:

            ```
            python PositionalError.py
            ```

    
    Round:
        Different Radius: simulations for different radius
        
        ## Plot mapping for apical length with respect to cell area and apical tension
        ```
            python ALMap.py
        ```
        ## Plot theoretical model for different radius showing apical length with respect to cell area
        ```
            python ALderived.py
        ```
        ## Plot final configurations
        ```
            python ConfigVisual.py
            python Plots.py
        ```
 ```           

## Run Lutolf simulation for Different Density (deterministic)

    ```
    sbatch --array=0-14 batchLutolfDensity
    ```
## Run Lutolf simulation for Different Density (noise)

    ```
    sbatch --array=0-4 batchLutolfNoiseDensity
    ```

## Run Lutolf simulation for Different Length (deterministic)

    ```
    sbatch --array=0-14 batchLutolfLength
    ```
## Run Lutolf simulation for Different Length (noise)

    ```
    sbatch --array=0-4 batchLutolfNoiseLength
    ```

## Run Lutolf simulation for Different Radius (deterministic)

    ```
    sbatch --array=0-14 batchLutolfRadius
    ```
## Run Lutolf simulation for Different Length (noise)

    ```
    sbatch --array=0-4 batchLutolfNoiseRadius
    ```

## Run Round simulation for Different Radius (deterministic)

    ```
    sbatch --array=0-4 batchLutolfNoiseRadius
    ```


```
Wavy
    ApicalTension: simulations for different values of apical tension
    SubstrateAmplitude: simulations for different values of substrate amplitude
    Wavelenght: simulations for different values of wavelenght               
```

## Run simulation

    ```
    python SimulateSinusoidal1D.py
    ```

    ## Plots:

    ```
    python Plots.py
    ```
