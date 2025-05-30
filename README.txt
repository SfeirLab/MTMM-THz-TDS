To run the code, simply open and run the main.m file.

Trying Different Structures:
To test different sample structures, modify the following variables in main.m using the data provided in the Data folder. 
1. SampleName – Specifies the sample for refractive index calculation. This should be a JSON file that defines the settings and structure of the sample.
2. Reference – Loads the experimental time trace of the reference measurement.
3. Sample – Loads the experimental time trace of the sample measurement.

Optional Settings:
For higher accuracy (with increased runtime) or faster execution (with potentially lower accuracy), you can adjust the following parameters:
PopSize – Population size. Increasing this improves accuracy but also increases runtime.
Maxit – Maximum number of generations. More iterations yield better results but take longer to compute.

You can also add your own data to the Data folder if needed. To do this:
1. Create a JSON file similar to the existing ones. You can add or remove layers by editing the thickness and refractive index values as needed.
2. Include two time-domain data files (as .txt files):
   - One for the reference measurement
   - One for the sample measurement
These files should contain the time-dependent experimental data and must also be placed in the Data folder.

This work was supported by the Gordon and Betty Moore Foundation, grant 10.37807/gbmf12235.

