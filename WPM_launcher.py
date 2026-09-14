import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import glob
from datetime import datetime
from helpers import make_gaussian, make_lens_slice, wpm_step, calc_sag, calc_overlap_1d

class LensGUI:
    def __init__(self, root):
        self.root = root
        root.title("Lens Sim WPM")
        
        frm_ctrl = tk.Frame(root)
        frm_ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        
        self.vars = {}
        # Map UI labels to internal dictionary keys for backward compatibility
        params = [
            ("Wavelength (μm)", "wl", 1.55), 
            ("Input beam w0 (μm)", "w0", 2.0), 
            ("R1 (μm)", "R1", 0.0), 
            ("k1", "k1", 0.0),
            ("R2 (μm)", "R2", -50.0), 
            ("k2", "k2", 0), 
            ("Lens center thickness (μm)", "thick", 50.0), 
            ("Lens diameter (μm)", "diam", 50.0),
            ("n_lens", "n_lens", 1.53), 
            ("n_background", "n_bg", 1.0), 
            ("Focal distance (μm)", "f_dist", 50.0), 
            ("Output beam w0 (μm)", "w_out", 5.2),
            ("Resolution", "res_div", 20.0) 
        ]
        
        for r, (label, name, val) in enumerate(params):
            tk.Label(frm_ctrl, text=label).grid(row=r, column=0, sticky='e')
            v = tk.DoubleVar(value=val)
            tk.Entry(frm_ctrl, textvariable=v, width=10).grid(row=r, column=1)
            self.vars[name] = v
            
        r += 1
        self.run_sim_var = tk.BooleanVar(value=True)
        tk.Checkbutton(frm_ctrl, text="Run Sim (Uncheck = Load)", variable=self.run_sim_var).grid(row=r, column=0, columnspan=2)
        
        r += 1
        tk.Button(frm_ctrl, text="GO", command=self.do_action, bg="green", fg="white", font=("Arial", 12, "bold")).grid(row=r, column=0, columnspan=2, pady=15)
        
        self.fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(1, 3, figsize=(14, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    def find_match(self):
        files = glob.glob("sim_*.npz")
        for f in files:
            try:
                data = np.load(f)
                match = True
                for key, var in self.vars.items():
                    if key in data:
                        if not np.isclose(data[key], var.get()):
                            match = False
                            break
                    else:
                        match = False
                        break
                if match:
                    return f
            except Exception:
                pass
        return None

    def do_action(self):
        match_file = self.find_match()
        
        if self.run_sim_var.get():
            if match_file:
                messagebox.showinfo("Found", "Sim already done. Loading now.")
                self.load_file(match_file)
            else:
                self.run_sim()
        else:
            if match_file:
                self.load_file(match_file)
            else:
                messagebox.showwarning("No Sim", "No file match. Run sim first.")

    def run_sim(self):
        wl = self.vars["wl"].get()
        w0 = self.vars["w0"].get()
        R1 = self.vars["R1"].get()
        k1 = self.vars["k1"].get()
        R2 = self.vars["R2"].get()
        k2 = self.vars["k2"].get()
        thick = self.vars["thick"].get()
        diam = self.vars["diam"].get()
        n_lens = self.vars["n_lens"].get()
        n_bg = self.vars["n_bg"].get()
        f_dist = self.vars["f_dist"].get()
        w_out = self.vars["w_out"].get()
        res_div = self.vars["res_div"].get()

        x_span = diam * 2.0
        # Truncate domain to user targeted focal limit
        z_span = thick + f_dist + 10.0 
        
        dx = dz = wl / res_div
        
        x = np.arange(-x_span/2, x_span/2, dx)
        z = np.arange(-5.0, z_span, dz)

        k0 = 2.0 * np.pi / wl
        kx = 2.0 * np.pi * np.fft.fftfreq(len(x), d=dx)
        window = np.exp(-((x/(0.9*x_span/2))**10))

        E = make_gaussian(x, w0)
        E_xz = []

        print("Processing wave propagation...")
        for z_cur in z:
            n_slice = make_lens_slice(x, z_cur, R1, k1, R2, k2, thick, diam, n_lens, n_bg)
            E = wpm_step(E, kx, k0, n_slice, dz, window)
            E_xz.append(E.copy())

        print("Simulation complete.")
        E_xz = np.array(E_xz)
        
        now = datetime.now().strftime("%Y%m%d_%H%M")
        save_name = f"sim_{now}.npz"
        
        save_dict = {'x': x, 'z': z, 'E_xz': E_xz}
        for key, var in self.vars.items():
            save_dict[key] = var.get()
            
        np.savez(save_name, **save_dict)
        print(f"Saved data to: {save_name}")
        
        self.draw_plots(z, x, E_xz, diam, R1, k1, R2, k2, thick, w_out)

    def load_file(self, filename):
        print(f"Loading: {filename}")
        data = np.load(filename)
        self.draw_plots(data['z'], data['x'], data['E_xz'], data['diam'], 
                        data['R1'], data['k1'], data['R2'], data['k2'], 
                        data['thick'], data['w_out'])

    def draw_plots(self, z, x, E_xz, diam, R1, k1, R2, k2, thick, w_out):
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()

        # Spatial shift to set lens exit face at z=0
        z_shifted = z - thick

        I_xz = np.abs(E_xz)**2
        self.ax1.imshow(I_xz.T, extent=[z_shifted[0], z_shifted[-1], x[0], x[-1]], aspect='auto', cmap='hot', origin='lower')
        self.ax1.set_title("Transverse Intensity (XZ Plane)")
        self.ax1.set_xlabel("z relative to lens exit (um)")
        self.ax1.set_ylabel("x (um)")

        lens_x = np.linspace(-diam/2, diam/2, 100)
        lens_z_front = calc_sag(np.abs(lens_x), R1, k1) - thick
        lens_z_back = calc_sag(np.abs(lens_x), R2, k2) 
        
        self.ax1.fill_betweenx(lens_x, lens_z_front, lens_z_back, color='cyan', alpha=0.3)
        self.ax1.plot([lens_z_front[0], lens_z_back[0]], [lens_x[0], lens_x[0]], color='cyan', linewidth=1.5)
        self.ax1.plot([lens_z_front[-1], lens_z_back[-1]], [lens_x[-1], lens_x[-1]], color='cyan', linewidth=1.5)

        eta_z = np.zeros(len(z))
        mfd_z = np.zeros(len(z))
        center_idx = len(x) // 2

        for i in range(len(z)):
            eta_z[i] = calc_overlap_1d(E_xz[i], x, w_out)
            
            slice_I = np.abs(E_xz[i])**2
            peak = slice_I[center_idx]
            if peak > 0:
                thresh = peak * np.exp(-2)
                idx_r = center_idx
                while idx_r < len(x)-1 and slice_I[idx_r] > thresh: 
                    idx_r += 1
                idx_l = center_idx
                while idx_l > 0 and slice_I[idx_l] > thresh: 
                    idx_l -= 1
                mfd_z[i] = np.abs(x[idx_r] - x[idx_l])    
        
        # Causality filter: true focusing only happens in free space post-lens (z_shifted > 0)
        valid_indices = np.where(z_shifted > 0)[0]
        if len(valid_indices) > 0:
            eta_z_valid = eta_z[valid_indices]
            max_idx = valid_indices[np.argmax(eta_z_valid)]
        else:
            max_idx = np.argmax(eta_z) # Failsafe
            
        max_eta = eta_z[max_idx]
        max_z_shifted = z_shifted[max_idx]
        optimal_mfd = mfd_z[max_idx]
        
        self.ax2.plot(z_shifted, eta_z, color='blue', linewidth=2)
        self.ax2.axvline(max_z_shifted, color='red', linestyle='--', label=f'Max Eta: {max_eta:.3f}\nf = {max_z_shifted:.2f} um')
        self.ax2.set_title(f"Coupling Efficiency (w_out = {w_out})")
        self.ax2.set_xlabel("z (um)")
        self.ax2.set_ylabel("Efficiency")
        self.ax2.legend()

        self.ax3.plot(z_shifted, mfd_z, color='purple', linewidth=2)
        self.ax3.axvline(max_z_shifted, color='red', linestyle='--', label=f'MFD at focal point: {optimal_mfd:.2f} um')
        self.ax3.set_title("MFD over Z")
        self.ax3.set_xlabel("z (um)")
        self.ax3.set_ylabel("MFD (um)")
        self.ax3.legend()

        self.fig.tight_layout()
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = LensGUI(root)
    root.mainloop()