# WPM_lens_simulator
Using WPM allows us for fast and accurate fiber/PIC lenses development

**USAGE GUIDE**

Parameter Configuration: Enter the desired optical and geometric parameters in the left control panel.

Run Sim / Load:

Keep "Run Sim" checked to compute a new propagation profile. If an exact match for your parameters already exists in the working directory (saved as a .npz file), the application will automatically load it to save time.

Uncheck "Run Sim" to force the application to look for an existing save file without computing.

Execution: Click GO. The simulation will execute and update the three plots on the right.

**Parameters**

Wavelength: Operating wavelength of the input beam (in $\mu$m).

Input beam w0: The beam radius ($1/e^2$) of the input Gaussian field before the lens.

R1 & k1: Radius of curvature and conic constant of the front lens surface.

R2 & k2: Radius of curvature and conic constant of the back lens surface.

Lens center thickness: On-axis thickness of the lens (in $\mu$m).

Lens diameter: Total physical aperture of the lens.

n_lens & n_background: Refractive indices of the lens material and the surrounding medium.

Focal distance: The expected focal distance. 

The simulation limits the spatial calculation to Lens center thickness + Focal distance + 10.0 $\mu$m to optimize compute time.

Output beam w0: The target beam radius ($1/e^2$) used to calculate the coupling efficiency (overlap integral).

Resolution: Grid resolution factor. Higher values yield finer spatial grids (dx = Wavelength / Resolution) but increase computation time.

**Assumptions & Physics**

Engine2D Dimensionality Reduction: To optimize computational load, the simulation evaluates a 1D transverse scalar field propagating along the Z-axis (an XZ slice). 

It assumes the beam and lens lack astigmatism and can be treated via 2D Cartesian coordinates.

Coupling Efficiency Approximation: Because the numerical engine computes a 1D transverse profile, the standard 1D overlap integral is squared during the final evaluation. This mathematically approximates the true 2D coupling efficiency for standard circularly symmetric optical fibers.

Z-Axis Translation: The coordinate system is dynamically shifted so that the back vertex of the lens sits exactly at $z = 0$. The spatial Z-coordinates in the output plots represent the free-space propagation distance directly.

Absorbing Boundaries: An exponential window function is applied at the edges of the computational grid to suppress non-physical reflections and aliasing from the Fast Fourier Transform (FFT) periodic boundaries.

**Known Limitations**

Focal Search Constraints: The algorithm is hard-coded to ignore any "virtual" focal points or localized intensity peaks inside the lens material. 

The maximum coupling efficiency and ideal MFD are strictly evaluated in the free-space region ($z > 0$).

Grid Spillage: If the simulated beam diverges significantly and hits the edges of the computational grid, the absorbing boundary window will artificially truncate the beam, which can skew the MFD and coupling efficiency calculations. 

Increase the lens diameter or decrease the propagation distance if the beam clips the vertical axis limits.

Scalar Approximation: The WPM engine relies on scalar diffraction theory. Vectorial effects, polarization dependencies, and significant anti-reflection (AR) coating behaviors are not modeled.

**Output Interpretation**

Transverse Intensity (XZ Plane): A heatmap of the beam intensity over space. The physical boundaries of the lens are overlaid in cyan. The zero-point on the X-axis represents the optical axis, and the zero-point on the Z-axis represents the exit face of the lens.

Coupling Efficiency: A plot of the overlap integral between the simulated beam and the target Output beam w0 across the Z-axis. The red dashed line pinpoints the exact Z-coordinate (focal length) that yields the maximum efficiency.

MFD over Z: The evolution of the Mode Field Diameter across the propagation axis, calculated dynamically by finding the $1/e^2$ intensity thresholds. The legend displays the specific MFD achieved at the focal plane determined by the efficiency plot.
