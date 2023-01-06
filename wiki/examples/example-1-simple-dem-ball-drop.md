# Example 1: Simple DEM Ball Drop

Follow the below instructions to compile and run a basic Project Chrono project on Quest. These instructions assume that you have a Northwestern **NetID** and access to the group project allocation **p31861**.&#x20;

This project will build and then run a wave tank of particles.

<figure><img src="../../.gitbook/assets/particles.JPG" alt=""><figcaption><p>Snapshot from Example 1 rendering in Paraview</p></figcaption></figure>

<details>

<summary>Step 1: Create a Working Directory</summary>

Create a directory for all of your developments and testing. No files/folders should be created or changed at the top-most '/projects/p31861' directory.

* Within FileZilla, enter your user directory folder&#x20;
* Right click in the folder and select **Create Directory and Enter It**
* Name the folder **workdir**

</details>

<details>

<summary>Step 2: Copy the Example Project </summary>

Copy the example project directory and files to your working directory

* Copy example project directory to your working directory, being sure to replace **LastnameFirstname** with your correct directory name&#x20;

```
cp -R /projects/p31861/ExampleProjects/exampleProject1 /projects/p31861/Users/LastnameFirstname/workdir 
```

* Navigate to the newly copied 'exampleProject1Make.sh' and 'exampleProject1Submit.sh' files in FileZilla and double-click on each to edit. Change all instances of **LastnameFirstname** in the files to your appropriate directory name and save/upload edited file back to Quest

</details>

<details>

<summary>Step 3: Compile and Run Project</summary>

In your SSH client navigate to the **exampleProject1** directory and run the following command to submit make job

```
sbatch exampleProject1Make.sh
```

When the job is complete, a build directory should have been created and inside will be an executable for your project. To run the executable, run the following command to submit run job

```
sbatch exampleProject1Submit.sh
```

_**Helpful Hint:** Reminder that the command 'squeue -u NetID' (with NetID swapped for your NetID) can be used to check the status of a job._

</details>
