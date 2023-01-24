# 💻 Desktop Installation and Setup

Follow the below instructions to get set up with the Chrono Preprocessor. These instructions were written with Windows users in mind. It is also possible to install on Mac or Linux.

Please note the versions of each piece of software. Newer or alternate versions may work, but have not been tested and verified for compatibility.

<details>

<summary>Step 1: Install FreeCAD</summary>

Install the latest version of FreeCAD (use at least version 0.20.2). The download is available for free:

[https://www.freecadweb.org/downloads.php](https://www.freecadweb.org/downloads.php)

</details>

<details>

<summary>Step 2: Install Qt Designer</summary>

Qt Designer normally ships as a part of Qt Creator. This is Qt's official editor and lets you do a lot more than just graphically design user interfaces. It is a full-fledged and very powerful C++ IDE. This power comes at a price however: The download for Qt Creator is gigabytes in size. You can download the full [Qt Creator suite](https://www.qt.io/product/development-tools). Otherwise, you can find alternative sites which provide installations with only Qt Designer. One such option is below:

[https://build-system.fman.io/qt-designer-download](https://build-system.fman.io/qt-designer-download)

</details>

<details>

<summary>Step 3: Install Git Client</summary>

Any Git client can be used to push and pull from the GitHub. We recommedn using SourceTree and these instructions will assume you are using that installation. You can download SourceTree for free here:

[https://www.sourcetreeapp.com/](https://www.sourcetreeapp.com/)

</details>

<details>

<summary>Step 4: Pull GitHub</summary>

We recommend pulling the GitHub directly into the FreeCAD workbench directory. Otherwise if you pull to another location then you will need to copy the pulled files to the appropriate directory.

* Open Sourcetree
* Select **File** > **Clone / New...**
* Select "**Remote**" and "**Add an account...**"
* For "Hosting Service" select "**GitHub**". For "Authentication" select "**OAuth**"
* Click on "**Refresh OAuth Token**" and login to GitHub and allow Sourcetree in the browser window that opens
* Click **Ok** in Sourcetree. Then "**chrono-preprocessor**" should populate on the right side of the window. If it doesn't, you may need to click refresh.
* Select "**chrono-preprocessor**" and click "**Clone**"
* When filling out the clone window, it should look like the image below, but with the appropriate username filled out instead of "**\<usr>**"
* Click "**Clone**"

</details>

<figure><img src="../../.gitbook/assets/clone.png" alt=""><figcaption></figcaption></figure>

<details>

<summary>Step 5: Test Installation</summary>

Verify that everything is installed properly by opening FreeCAD and check if the Chrono Workbench is available in the list of installed workbenches.

</details>
