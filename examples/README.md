# Energy_calculation
<h2>StackGen example problems</h2>

There are 14 sub-directories present here. Each sub-directory contains a sample problem you can run with StackGen. 
Each problem has a user input script (input.user), DFTB+ input file (dftb_in.hsd), PSO trajectories files (swarm_traj.out, traj.out)
as output, optimized dimer structure (generated_supramo_2.xyz) and three image files of monomer configuration, side view, and top view of 
energy-optimized dimer structure, respectively.

The following molecules are considered sample problems. The input and output files are present in their corresponding sub-directories:
<ol>
  <li>BTA_68b</li>
  <li>BTA_66d </li>
  <li>Triazine-trisamide-1 (TTA-1) </li>
  <li>Triazine-trisamide-2 (TTA-2)</li>
  <li>PDI_20</li>
  <li>PDI_92</li>
  <li>PDI_branch</li>
  <li>PDI_vshape</li>
  <li>Hexabenzocoronene(HBC)</li>
  <li>BTA</li>
  <li>Ph-PDI</li>
  <li>Alkyl-PDI</li>
  <li>PDI-1-C1</li>
  <li>PBI-1 </li>
   
</ol>  

Here's an example of how to run and visualize the energy-optimized structures of these sample problems:


cd BTA_66d </br>
source ../../path_to_set_initenv.sh file </br>
stackgen_run  </br>

Running the simulation produces two PSO trajectory files, namely traj.out and swarm_traj.out. traj.out captures the progress of the best fitness value found by the entire swarm at each iteration. Each line in traj.out corresponds to an iteration, detailing the best fitness value and the order parameter values (tx, ty, tz, twist) at that iteration. Whereas swarm traj.out records the complete trajectory of the PSO process. It documents the fitness value and the position of swarm particles for each iteration.
The energy-optimized dimer structure can be visualized using VMD.

<i>StackGen</i> generated order parameters of the energy-minimized configuration of these sample problems are depicted below:
<table>
 <tr>
    <th>Molecule</th>
    <th>Dimension</th>
    <th>tx</th>
    <th>ty</th>
    <th>tz</th>
    <th>&#952;</th>
    <tr>
    <td>BTA_68b</td>
    <td>4</td>
    <td>0.001</td>
    <td>0.332	</td>
    <td>3.478</td>
    <td>-65.428</td>
  </tr>
   <tr>
    <td>BTA_66d</td>
    <td>4</td>
    <td>0.552	</td>
    <td>0.310</td>
    <td>3.588</td>
    <td>-57.817</td>
  </tr>
  <tr>
    <td>Triazine-trisamide-1 (TTA-1)</td>
    <td>4</td>
    <td>0.184</td>
    <td>0.044</td>
    <td>3.291/td>
    <td>-31.013<</td>
  </tr>
  <tr>
    <td>Triazine-trisamide-2 (TTA-2)</td>
    <td>4</td>
    <td>0.287</td>
    <td>0.000 </td>
    <td>3.218</td>
    <td>-31.295</td>
  </tr>
   <tr>
    <td>PDI_20</td>
    <td>4</td>
    <td> 3.491	</td>
    <td> 0.357</td>
    <td>3.580</td>
    <td>0.003</td>
  </tr>
  <tr>
    <td>PDI_92</td>
    <td>4</td>
    <td>1.968	</td>
    <td> 0.580</td>
    <td>3.343</td>
    <td>-37.069</td>
  </tr>
   <tr>
    <td>PDI_branch</td>
    <td>4</td>
    <td> 1.505</td>
    <td>0.909 </td>
    <td>3.304</td>
    <td>-32.655</td>
  </tr>
    <tr>
    <td>PDI_vshape</td>
    <td>4</td>
    <td>1.035	</td>
    <td> 0.520</td>
    <td>3.438</td>
    <td>-39.749</td>
  </tr>
   <tr>
    <td>Hexabenzocoronene(HBC)</td>
    <td>4</td>
    <td>1.341	</td>
    <td>0.545</td>
    <td>3.305	</td>
    <td>-59.200</td>
  </tr>
   <tr>
    <td>BTA</td>
    <td>2</td>
    <td>0.0 </td>
    <td>0.0 </td>
    <td>3.427</td>
    <td>296.642</td>
  </tr>
   <tr>
    <td>Ph-PDI</td>
    <td>2</td>
    <td>0.0</td>
    <td>0.0</td>
    <td>3.388</td>
    <td>349.073</td>
  </tr>
    <tr>
    <td>Alkyl-PDI </td>
    <td>2</td>
    <td>0.0</td>
    <td>0.0 </td>
    <td>3.487  </td>
    <td>29.872</td>
  </tr>
  <tr>
    <td>PDI-1-C1</td>
    <td>4</td>
    <td>1.302 </td>
    <td>0.218</td>
    <td>3.314</td>
    <td>-20.991</td>
  </tr>
  <tr>
    <td>PBI-1</td>
    <td>4</td>
    <td>2.165</td>
    <td>0.030</td>
    <td>3.401</td>
    <td>19.988</td>
  </tr>

   
</tr>
</table>
