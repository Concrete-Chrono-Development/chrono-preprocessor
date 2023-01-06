# 💻 Desktop Installation and Setup

Follow the below instructions to get set up with Project Chrono. These instructions were written with Windows users in mind. It is also possible to install on Mac or Linux; Mac instructions are located in the notes section of each Step.

Please note the versions of each piece of software. Newer or alternate versions may work, but have not been tested and verified for compatibility.

<details>

<summary>Step 1: Install C++ Compiler</summary>

Install Microsoft Visual Studio 2022. The [community edition of the latest Visual Studio](https://visualstudio.microsoft.com/downloads/) is available for free.

* During installation make sure to check and include _"Desktop Development with C++"_
* After installation open Visual Studio and sign-in if necessary

```
- For Mac: Use Xcode Package. Download via App Store for free - it contains the clang++ compiler.
- Notes: Other compilers were also tested (e.g. Intel C++, PGI) but they are not officially supported and maintained. While it is likely possible to build Chrono with other toolchains, this might require changes to the CMake scripts.
- Notes: Any version of Visual Studio 2019 or newer should work. Visual Studio 2017 has problems with the heavy use of inlining in recent version of Eigen. For the latest version of Chrono (specifically due to the reimplemented ANCF elements), this can result in very long compilation times or even hang ups. We recommend using VS 2019 or newer.
```

</details>

<details>

<summary>Step 2: Install Eigen Library</summary>

Download the Eigen version 3.4.0 zipped  source code. This [code is available for free](https://gitlab.com/libeigen/eigen/-/releases/3.4.0).&#x20;

* Unzip the downloaded file and store the contents in an easy-to-find location, suggested location is: C:\workspace\libraries\eigen-3.4.0

<pre><code><strong>- For Mac: Install it via homebrew: brew install eigen. Homebrew installs into /opt/homebrew since MacOS 12 Monterey and the new Apple Silicon (arm46, M1, M2...) hardware. If Eigen is not found automatically, you can search its folder with:
</strong>find /opt/homebrew -name Eigen
<strong>- Notes: Any version of Eigen 3.3.0 or newer should work. 
</strong></code></pre>

</details>

<details>

<summary>Step 3: Install Irrlicht Library</summary>

Download Irrlicht SDK version 1.8.5. This [code is available for free](https://irrlicht.sourceforge.io/?page\_id=10).&#x20;

* Unzip the downloaded file and store the contents in an easy-to-find location, suggested location is: C:\workspace\libraries\irrlicht-1.8.5

<pre><code><strong>- For Mac: The best way to install irrlicht on the Mac is: brew install irrlicht (release v.1.8.5). On MacOS 12 (Monterey) you have to set IRRLICHT_ROOT to /opt/homebrew.
</strong><strong>- Notes: Any version of Irrlicht SDK version 1.8.2 or newer should work. 
</strong></code></pre>

</details>

<details>

<summary>Step 4: Install CMake</summary>

Install CMake version 3.25.0. An installer for the [software is available for free](https://cmake.org/download/).

* During installation be sure to check for the CMake executable to be included in your Path environmental variable

<pre><code><strong>- For Mac: The CMake.app bundle also contains command line tools, you must set appropriate links to use it from the terminal. It is better to install a pure command line version via homebrew (https://brew.sh). After installing the home brew package manager type: brew install cmake in the terminal.
</strong><strong>- Notes: Any version of CMake version 1.8.2 or newer should work. 
</strong></code></pre>

</details>

<details>

<summary>Step 5: Install Git Client - In Progress</summary>

*

</details>

<details>

<summary>Step 6: Clone Project Chrono - <em>In Progress</em></summary>

*

</details>

<details>

<summary>Step 7: Run CMake - <em>In Progress</em></summary>

*

</details>

<details>

<summary>Step 8: Compile Project - <em>In Progress</em></summary>

*

</details>

Source for many of these instructions is Project Chrono documentation. More details [here.](https://api.projectchrono.org/tutorial\_install\_chrono.html)
